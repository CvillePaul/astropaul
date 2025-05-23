"""
SkyNet API: object-relational mapping

This module constitutes the object-oriented API interface. It contains a set of
classes that are counterparts to SQLAlchemy ORM classes defined in
skynet.db.model and can be used in the same way on the client side or within
SkyNet web2py controllers. API ORM objects may be treated as "static" versions
of SQLAlchemy ORM objects that do not require a database connection and are thus
more efficient; whenever the actual database access is required, controller
calls the consolidated code (see :mod:`skynet.api.handlers`), which gets API ORM
objects on input to update the database and/or returns them on output as a
result of a database query.

Serialization to/deserialization from dictionaries, JSON strings, or
:mod:`skynet.db.model` objects is done via marshmallow schemas in
:mod:`skynet.api.serializers`. There is a one-to-one mapping between
:mod:`skynet.db.model` classes, :mod:`skynet.api.serializers` schemas, and
:mod:`skynet.api.orm` objects.
"""

import math
from datetime import datetime, timedelta
import pprint

from .serializers import *
from .client import AttrDict
from .methods import (
    collabs as collab_methods, exps as exp_methods, filters as filter_methods,
    groups as group_methods, jobs as job_methods,
    master_cals as master_cal_methods, obs as obs_methods,
    radio_obs as radio_obs_methods, scopes as tele_methods,
    triggers as trigger_methods, users as user_methods)
from .errors import BadRequestError, InternalApiError


__all__ = [
    'Antenna',
    'Camera',
    'CampaignTrigger',
    'Collaboration',
    'DarkCalibrationConfig',
    'DataProcessingFile',
    'DataProcessingInputFile',
    'DataProcessingJob',
    'DataProcessingOutputFile',
    'DmsStatus',
    'Exposure',
    'Filter',
    'FilterTransmission',
    'FlatCalibrationConfig',
    'GenericFilter',
    'Group',
    'Job',
    'LocalHorizon',
    'MasterCalibration',
    'MasterComponent',
    'Observation',
    'ObsNote',
    'ObsStats',
    'ObsUserPreference',
    'OrbitalElements',
    'PendingCollabMembership',
    'PendingGroupMembership',
    'PriorityAccessSchedule',
    'RadioCartographerJob',
    'RadioFilter',
    'RadioObs',
    'RadioObsDaisy',
    'RadioObsMap',
    'RadioObsOnOff',
    'RadioObsTrack',
    'Receiver',
    'ReceiverHistory',
    'Role',
    'Site',
    'StandardFilter',
    'TeleOwner',
    'TeleOwnerShare',
    'Telescope',
    'TelescopeBase',
    'TeleSiteRole',
    'TeleStatus2',
    'TimeAccount',
    'TimeAccountSettings',
    'TimeAccountTxn',
    'User',
    'UserSelfRegistration',
    'WcsSolution',
]


pp = pprint.PrettyPrinter(indent=4)


def relationship(field_name, orm_class, uselist=True):
    """
    Define a relationship property for an ORM class

    class Group(SkynetApiBase):
        ...
        owningCollab = relationship('owningCollabId', Collaboration)

    :param str field_name: name of the field corresponding to the relationship
    :param orm_class: API ORM class or class name to return
    :param bool uselist: list-type property (default) vs scalar property

    :return: new property
    """
    def getter(self):
        self.refresh()

        if isinstance(orm_class, str):
            cls = globals()[orm_class]
        else:
            cls = orm_class

        if cls.__get_method__ is None:
            # .get() is not supported by relationship class, return the original
            # field value as is; obtain field name by comparing the class'
            # property getters to this function
            for name in dir(self.__class__):
                val = getattr(self.__class__, name)
                if isinstance(val, property) and val.fget is getter:
                    val = self.fields[name]
                    if isinstance(val, dict):
                        val = AttrDict(val)
                    return val
            return None

        try:
            val = self.fields[field_name]
        except KeyError:
            return None
        if val is None:
            return None

        if uselist:
            return [cls.get(item, api_key=self._api_key) for item in val]
        return cls.get(val, api_key=self._api_key)

    return property(getter)


class SkynetApiBase(object):
    """
    Base class for SkyNet API ORM objects

    ORM objects are defined as follows:
        from skynet.api.serializers import MySchema
        from skynet.api.methods import my_methods

        class MyClass(SkynetApiBase):
            __schema__ = MySchema
            __get_method__ = staticmethod(my_methods.get)
            __query_method__ = staticmethod(my_methods.query)
            __update_method__ = staticmethod(my_methods.update)

            # Any custom methods and properties follow; you may use them to
            # mimic the original skynet.db.model class behavior

    Virtual constructors (i.e. returning different subclasses of the given class
    based on the value of a certain field supplied on creation) are supported as
    follows:

        class MyVirtualClass(SkynetApiBase):
            ...
            __polymorphic_on__ = 'field'
            __polymorphic_table__ = {'value1': 'Subclass1',
                                     'value2': 'Subclass2',
                                     ...
                                    }

        class Sublass1(MyVirtualClass):
            ...

        class Sublass2(MyVirtualClass):
            ...

    Usage:
        obj = MyClass(<SQLAlchemy ORM class instance>)
        obj = MyClass('JSON-serialized skynet.api.serializers schema')
        obj = MyClass({'field': value, 'field': value...})
        obj = MyClass(field=value, field=value...)
        obj = MyClass.get(id_or_name)
        objects = MyClass.query(param=value...)

        # Create an instance of the corresponding subclass based on the
        # polymorphic field value (see e.g. RadioObs); for non-polymorphic
        # classes, MyClass.create(...) is the same as MyClass(...)
        obj = MyVirtualClass.create(<SQLAlchemy ORM class instance>)
        obj = MyVirtualClass.create('JSON-serialized schema')
        obj = MyVirtualClass.create({'field': value, 'field': value...})
        obj = MyVirtualClass.create(field=value, field=value...)

        print(obj.field) or print(obj['field'])  # get the current value
        obj.field = ... or obj['field'] =  # assign a new value with validation
        del obj.field or del obj['field']  # mark the field unspecified

        obj.update()  # commit changes to object fields to the database
        obj.refresh()  # update static fields from the database

        obj.dump()  # serialize object to dictionary
        obj.dump(skynet.db.model.MyClass)  # serialize to SQLAlchemy ORM object
        obj.dumps()  # serialize to JSON string

    To authenticate the possible client-side procedural API calls made by
    custom methods and properties, the user may provide an API key on creation
    of an ORM object. This is only required when the object is instantiated by
    a website controller. For the client side, the key is set automatically.
    """
    __schema__ = None  # marshmallow schema class
    __get_method__ = None  # skynet.api.orm.methods.*.get() method reference
    __query_method__ = None  # skynet.api.orm.methods.*.query() method reference
    __update_method__ = None  # skynet.api.orm.methods.*.update()
    __polymorphic_on__ = None
    __polymorphic_table__ = None
    schema = None  # schema class instance
    fields = None  # dictionary of field values
    _api_key = None  # API key used to authenticate the client API calls

    def __init__(self, _object=None, _api_key=None, **kwargs):
        """
        Initialize an API ORM object from an object, a JSON string, or a set of
        field=value pairs

        :param dict | str | object _object: an object (e.g. a skynet.db.model
            class instance), a dictionary, or or a JSON string to initialize
            from; all fields are serialized using the object's marshmallow
            schema; if omitted, initialize from keyword=value pairs given as
            keyword arguments
        :param str _api_key: optional API access token; API clients may use this
            to authenticate as a different user than the one stored in
            ~/.skynet/API_KEY.
        :param kwargs: a set of field-value pairs; if both `_object` and
            field=value pairs are present, the latter override the former
        """
        schema_class = object.__getattribute__(self, '__schema__')
        if schema_class is None:
            raise InternalApiError(
                'Missing schema for SkyNet API ORM class "{}"'.format(
                    object.__getattribute__(self, '__class__').__name__))
        if not issubclass(schema_class, Schema):
            raise InternalApiError(
                '__schema__ for SkyNet API ORM class "{}" is not a '
                'marshmallow.Schema'.format(
                    object.__getattribute__(self, '__class__').__name__))
        schema = schema_class()
        object.__setattr__(self, 'schema', schema)
        object.__setattr__(self, '_api_key', _api_key)

        if _object is not None:
            if isinstance(_object, int):
                # Initializing from an integer ID; handle below
                _object = str(_object)

            if isinstance(_object, str) or isinstance(_object, type(u'')):
                # Load from a JSON string
                try:
                    # noinspection PyTypeChecker
                    _object = schema.loads(_object)
                except (ValidationError, ValueError) as e:
                    # noinspection PyUnresolvedReferences
                    if isinstance(e, ValidationError) and (
                            hasattr(e, 'field_names') and  # marshmallow 2
                            schema in e.field_names or
                            hasattr(e, 'field_name') and  # marshmallow 3
                            e.field_name == '_schema') or \
                            isinstance(e, ValueError):
                        # Initializing from object ID or name
                        try:
                            _object = int(_object)
                        except ValueError:
                            # Treat input object as ORM object name
                            try:
                                name_field = schema.__name_field__
                            except AttributeError:
                                name_field = 'name'
                            _object = {name_field: _object}
                        else:
                            # Treat input object as ORM object ID
                            try:
                                id_field = schema.__id_field__
                            except AttributeError:
                                id_field = 'id'
                            _object = {id_field: _object}
                    else:
                        raise

            elif isinstance(_object, SkynetApiBase):
                # When initializing from an API ORM object, get its fields
                # directly
                _object = object.__getattribute__(_object, 'fields')

            # Get all fields declared in the schema
            all_fields = schema.dump(_object)

            # Override with the explicit field=value pairs
            all_fields.update(schema.dump(kwargs))
        else:
            # Load from a dictionary of keyword=value pairs; validate all
            # keywords first
            all_fields = schema.dump(kwargs)

        # Make sure the polymorphic selector field has the appropriate value
        polymorphic_field = self.__polymorphic_on__
        if polymorphic_field:
            for name, val in self.__polymorphic_table__.items():
                cls = object.__getattribute__(self, '__class__')
                if isinstance(val, str):
                    try:
                        val = globals()[val]
                    except KeyError:
                        raise InternalApiError(
                            'Unknown subclass "{}" in "{}" polymorphic '
                            'table'.format(val, cls.__name__))
                if cls is val:
                    all_fields[polymorphic_field] = name
                    break

        # Initialize API ORM object fields
        object.__setattr__(
            self, 'fields', schema.load(all_fields))

    @classmethod
    def create(cls, _object=None, _api_key=None, **kwargs):
        """
        SkynetApiBase virtual constructor; creates an API object instance of the
        appropriate subtype based on the supplied polymorphic field value

        :param dict | str | object _object: an object, a dictionary, or or a
            JSON string to initialize from
        :param str _api_key: optional API access key
        :param kwargs: a set of field=value pairs

        :return: instance of the appropriate API class
        :rtype: SkynetApiBase
        """
        polymorphic_field = cls.__polymorphic_on__
        if not polymorphic_field:
            return cls(_object, _api_key, **kwargs)

        polymorphic_identity = None
        if _object is not None:
            if isinstance(_object, str) or isinstance(_object, type(u'')):
                try:
                    _object = json.loads(_object)
                except ValueError:
                    # Initializing from ID/name, polymorphic field undefined,
                    # instantiate the base class
                    return cls(_object, _api_key, **kwargs)
            elif isinstance(_object, SkynetApiBase):
                _object = object.__getattribute__(_object, 'fields')

            try:
                polymorphic_identity = getattr(_object, polymorphic_field)
            except AttributeError:
                try:
                    polymorphic_identity = _object[polymorphic_field]
                except (AttributeError, KeyError, TypeError):
                    pass
        try:
            polymorphic_identity = kwargs[polymorphic_field]
        except KeyError:
            pass

        try:
            subclass = cls.__polymorphic_table__[polymorphic_identity]
        except (KeyError, TypeError):
            subclass = cls
        else:
            if isinstance(subclass, str):
                try:
                    subclass = globals()[subclass]
                except KeyError:
                    raise InternalApiError(
                        'Unknown subclass "{}" in "{}" polymorphic '
                        'table'.format(subclass, cls.__name__))

        return subclass(_object, _api_key, **kwargs)

    def __getattribute__(self, item):
        """
        Get an API ORM object attribute or field value

        :param item: attribute/field name

        :return: attribute/field value
        :rtype: None | int | str | float | datetime | dict | list
        """
        try:
            return object.__getattribute__(self, item)
        except AttributeError:
            if item in object.__getattribute__(self, 'schema').fields:
                try:
                    return object.__getattribute__(self, 'fields')[item]
                except KeyError:
                    raise AttributeError(
                        '"{}" object has no field "{}"'.format(
                            object.__getattribute__(self, '__class__').__name__,
                            item))
            else:
                raise
    __getitem__ = __getattribute__

    def __setattr__(self, key, value):
        """
        Set attribute or field value; validates the value being assigned to
        known fields, raises AttributeError for non-field attributes

        :param key: attribute/field name
        :param value: new attribute/field value

        :rtype: None
        """
        schema = object.__getattribute__(self, 'schema')
        if key in schema.fields:
            field_vals = object.__getattribute__(self, 'fields')
            new_val = schema.load(schema.dump({key: value}))[key]

            # Don't allow changing the polymorphic selector field unless we are
            # a base class instance
            if key == object.__getattribute__(
                    self, '__polymorphic_on__') and key in field_vals and \
                    value != field_vals[key]:
                cls = object.__getattribute__(self, '__class__')
                for subclass in object.__getattribute__(
                        self, '__polymorphic_table__').values():
                    if isinstance(subclass, str):
                        try:
                            subclass = globals()[subclass]
                        except KeyError:
                            raise InternalApiError(
                                'Unknown subclass "{}" in "{}" polymorphic '
                                'table'.format(subclass, cls.__name__))
                    if cls is subclass:
                        raise AttributeError(
                            'Cannot change object "{}" field "{}"'.format(
                                cls.__name__, key))

            field_vals[key] = new_val

        else:
            raise AttributeError(
                'Cannot set object "{}" attribute "{}"'.format(
                    object.__getattribute__(self, '__class__').__name__, key))

    __setitem__ = __setattr__

    def __delattr__(self, item):
        """
        Delete attribute or field value; raises AttributeError for non-field
        attributes

        :param item: attribute/field name

        :rtype: None
        """
        if item in object.__getattribute__(self, 'schema').fields:
            # Don't allow deleting the polymorphic selector field unless we are
            # a base class instance
            if item == object.__getattribute__(self, '__polymorphic_on__'):
                cls = object.__getattribute__(self, '__class__')
                for subclass in object.__getattribute__(
                        self, '__polymorphic_table__').values():
                    if isinstance(subclass, str):
                        try:
                            subclass = globals()[subclass]
                        except KeyError:
                            raise InternalApiError(
                                'Unknown subclass "{}" in "{}" polymorphic '
                                'table'.format(subclass, cls.__name__))
                    if cls is subclass:
                        raise AttributeError(
                            'Cannot delete object "{}" field "{}"'.format(
                                cls.__name__, item))

            try:
                del object.__getattribute__(self, 'fields')[item]
            except KeyError:
                pass

        else:
            raise AttributeError(
                'Cannot delete object "{}" attribute "{}"'.format(
                    object.__getattribute__(self, '__class__').__name__, item))
    __delitem__ = __delattr__

    def __contains__(self, item):
        """
        Does an API ORM object have the given field set?

        :param item: field name

        :return: True is the given field is set
        :rtype: bool
        """
        return item in object.__getattribute__(self, 'fields')

    def __str__(self):
        """
        Return a string representation of the API ORM object

        :return: pretty-printed object fields
        :rtype: str
        """
        return pp.pformat(self.fields)

    def __repr__(self):
        """
        Return an internal representation of the API ORM object

        :return: repr(obj)
        :rtype: str
        """
        return '<{} at 0x{:x} {}'.format(
            self.__class__.__name__, id(self), repr(self.fields))

    def dump(self, _class=None):
        """
        Serialize an API ORM object to a dictionary of field=value pairs or
        to another object (e.g. skynet.db.model class)

        :param _class: if given, serialize to the given class; otherwise,
            serialize to a dictionary

        :return: an instance of _class if specified or a dictionary otherwise
        :rtype: _class | dict
        """
        all_fields = object.__getattribute__(self, 'fields')
        if _class is None:
            return object.__getattribute__(self, 'schema').dump(all_fields)

        try:
            return _class(**all_fields)
        except TypeError:
            obj = _class()
            for name, val in all_fields.items():
                setattr(obj, name, val)
            return obj

    def dumps(self):
        """
        Serialize an API ORM object to a JSON string

        :return: serialized API ORM object
        :rtype: str
        """
        return object.__getattribute__(self, 'schema').dumps(
            object.__getattribute__(self, 'fields'))

    @classmethod
    def get(cls, id_or_name, api_key=None):
        """
        Return an API ORM object from the corresponding get() call, a shortcut
        to APIClass(skynet.api.methods.*.get(id_or_name))

        :param str | int id_or_name: object ID or name
        :param str api_key: optional API access token

        :return: API ORM object
        :rtype: SkynetApiBase
        """
        get_method = cls.__get_method__
        if get_method is None:
            raise BadRequestError(
                'Class "{}" does not support get()'.format(cls.__name__))

        return cls.create(
            get_method(id_or_name, api_key=api_key), _api_key=api_key)

    @classmethod
    def query(cls, api_key=None, **kwargs):
        """
        Return a list of API ORM objects returned by the corresponding query()
        call, a shortcut to
        [APIClass(obj) for obj in skynet.api.methods.*.query(**kwargs))]

        :param str api_key: optional API access token
        :param kwargs: query arguments

        :return: list of API ORM objects
        :rtype: list[SkynetApiBase]
        """
        query_method = cls.__query_method__
        if query_method is None:
            raise BadRequestError(
                'Class "{}" does not support query()'.format(cls.__name__))

        kwargs = dict(kwargs)
        kwargs['include'] = '*'
        kwargs['exclude'] = ''
        return [cls.create(obj, api_key=api_key)
                for obj in query_method(api_key=api_key, **kwargs)]

    def refresh(self, _from=None):
        """
        Refresh static fields by making a skynet.api.methods.*.get() query or
        from a result of some other API call

        :param dict _from: result of an API call to set fields from (e.g.
            skynet.methods.user.create()); default: refresh from
            skynet.methods.*.get(self.id)
        :return: None
        """
        if _from is not None:
            object.__setattr__(self, 'fields', _from)

        get_method = object.__getattribute__(self, '__get_method__')
        if get_method is None:
            return

        try:
            obj_id = object.__getattribute__(self, 'fields')['id']
        except KeyError:
            # Object not yet initialized
            return

        res = get_method(
            obj_id, api_key=object.__getattribute__(self, '_api_key'))
        object.__setattr__(
            self, 'fields',
            object.__getattribute__(self, 'schema').load(res))

    def update(self):
        """
        Commit locally changed static object fields to the database by making
        a skynet.api.methods.*.update() query; this is an inverse to refresh()

        :return: None
        """
        update_method = object.__getattribute__(self, '__update_method__')
        if update_method is None:
            return

        all_fields = object.__getattribute__(self, 'schema').fields
        update_method(
            object.__getattribute__(self, 'fields')['id'],
            api_key=object.__getattribute__(self, '_api_key'),
            **{name: val
               for name, val in object.__getattribute__(self, 'fields').items()
               if is_scalar_field(all_fields[name])})


class CampaignTrigger(SkynetApiBase):
    __schema__ = CampaignTriggerSchema
    __get_method__ = staticmethod(trigger_methods.get)
    __query_method__ = staticmethod(trigger_methods.query)
    __update_method__ = staticmethod(trigger_methods.update)

    obs = relationship('obsId', 'Observation', uselist=False)
    refFilter = relationship('filterId', 'StandardFilter')
    refFilterRequested = relationship('refFilterId', 'Filter', uselist=False)
    campaign = relationship('campaignId', 'Campaign')

    def submit(self):
        """
        Submit a trigger with parameters set by the current instance fields

        :return: None; the trigger is updated in place
        :rtype: None
        """
        self.refresh(trigger_methods.add(api_key=self._api_key, **self.fields))

    def activate(self):
        trigger_methods.update(
            trigger=self.id, state='active', api_key=self._api_key)
        self.refresh()

    def cancel(self):
        trigger_methods.update(
            trigger=self.id, state='canceled', api_key=self._api_key)
        self.refresh()


class Camera(SkynetApiBase):
    __schema__ = CameraSchema


class Collaboration(SkynetApiBase):
    __schema__ = CollaborationSchema
    __get_method__ = staticmethod(collab_methods.get)
    __query_method__ = staticmethod(collab_methods.query)
    __update_method__ = staticmethod(collab_methods.update)

    groups = relationship('groups', 'Group')
    timeAccounts = relationship('timeAccounts', 'TimeAccount')

    @property
    def entityId(self):
        return 'c{:d}'.format(self.id)

    @property
    def entityType(self):
        return 'Collaboration'

    @property
    def ownedGroups(self):
        # noinspection PyTypeChecker
        return [g for g in self.groups
                if hasattr(g, 'owningCollabId') and g.owningCollabId == self.id]

    pendingCollabMemberships = relationship(
        'pendingCollabMemberships', 'PendingCollabMembership')
    admins = relationship('admins', 'User')

    # TODO: def invite(self, group):

    # TODO: def requestMembership(self, group):


class DarkCalibrationConfig(SkynetApiBase):
    __schema__ = DarkCalibrationConfigSchema


class DataProcessingFile(SkynetApiBase):
    __schema__ = DataProcessingFileSchema
    __polymorphic_on__ = 'type'
    __polymorphic_table__ = {
        'input': 'DataProcessingInputFile',
        'output': 'DataProcessingOutputFile',
    }

    job = relationship('job', 'DataProcessingJob', uselist=False)


class DataProcessingInputFile(DataProcessingFile):
    pass


class DataProcessingOutputFile(DataProcessingFile):
    pass


class DataProcessingJob(SkynetApiBase):
    __schema__ = DataProcessingJobSchema
    __polymorphic_on__ = 'type'
    __polymorphic_table__ = {
        'radio_cartographer': 'RadioCartographerJob',
    }

    files = relationship('files', 'DataProcessingFile')


class DmsStatus(SkynetApiBase):
    __schema__ = DmsStatusSchema


class Exposure(SkynetApiBase):
    __schema__ = ExposureSchema
    __get_method__ = staticmethod(exp_methods.get)
    __query_method__ = staticmethod(exp_methods.query)
    __update_method__ = staticmethod(exp_methods.update)

    targetExp = relationship('targetExpId', 'Exposure', uselist=False)
    linkedExps = relationship('linkedExps', 'Exposure')
    filterRequested = relationship('filterIdRequested', 'Filter', uselist=False)
    filterUsed = relationship('filterIdUsed', 'StandardFilter', uselist=False)
    telescope = relationship('teleId', 'Telescope', uselist=False)
    telescopeRequested = relationship(
        'teleIdRequested', 'Telescope', uselist=False)
    timeAccount = relationship('timeAccountId', 'TimeAccount', uselist=False)
    parentGroupTimeAccount = relationship(
        'parentGroupTimeAccountId', 'TimeAccount', uselist=False)
    parentCollabTimeAccount = relationship(
        'parentCollabTimeAccountId', 'TimeAccount', uselist=False)
    teleOwner = relationship('teleOwnerIdUsed', 'TeleOwner', uselist=False)
    wcsSolution = relationship('wcsId', 'WcsSolution', uselist=False)
    obs = relationship('obsId', 'Observation', uselist=False)
    camera = relationship('camera', 'Camera', uselist=False)

    def submit(self, obs):
        """
        Add the current exposure to an existing observation

        :param int | str | skynet.api.orm.Observation obs: ID or API ORM object
            for the observation to add the exposure to

        :rtype: None
        """
        if hasattr(obs, 'id'):
            obs = obs.id
        self.refresh(exp_methods.add(
            obs=obs, api_key=self._api_key, **self.fields))

    def activate(self):
        exp_methods.update(exp=self.id, state='ready', api_key=self._api_key)
        self.refresh()

    def cancel(self):
        exp_methods.update(exp=self.id, state='canceled', api_key=self._api_key)
        self.refresh()


class Filter(SkynetApiBase):
    __schema__ = FilterSchema
    __get_method__ = staticmethod(filter_methods.get)
    __query_method__ = staticmethod(filter_methods.query)
    __polymorphic_on__ = 'type'
    __polymorphic_table__ = {
        'standard': 'StandardFilter',
        'generic': 'GenericFilter',
    }


class FilterTransmission(SkynetApiBase):
    __schema__ = FilterTransmissionSchema


class FlatCalibrationConfig(SkynetApiBase):
    __schema__ = FlatCalibrationConfigSchema


class GenericFilter(Filter):
    __schema__ = GenericFilterSchema


class Group(SkynetApiBase):
    __schema__ = GroupSchema
    __get_method__ = staticmethod(group_methods.get)
    __query_method__ = staticmethod(group_methods.query)
    __update_method__ = staticmethod(group_methods.update)

    owningCollab = relationship(
        'owningCollabId', 'Collaboration', uselist=False)

    users = relationship('users', 'User')

    @property
    def ownedUsers(self):
        # noinspection PyTypeChecker
        return [u for u in self.users
                if hasattr(u, 'owningGroupId') and u.owningGroupId == self.id]

    timeAccounts = relationship('timeAccounts', 'TimeAccount')
    collabs = relationship('collabs', 'Collaboration')
    pendingGroupMemberships = relationship(
        'pendingGroupMemberships', 'PendingGroupMembership')
    pendingCollabMemberships = relationship(
        'pendingCollabMemberships', 'PendingCollabMembership')
    admins = relationship('admins', 'User')

    @property
    def entityId(self):
        return 'g{:d}'.format(self.id)

    @property
    def entityType(self):
        return 'Group'


class Job(SkynetApiBase):
    __schema__ = JobSchema
    __get_method__ = staticmethod(job_methods.get)
    __query_method__ = staticmethod(job_methods.query)

    user = relationship('userId', 'User', uselist=False)


class LocalHorizon(SkynetApiBase):
    __schema__ = LocalHorizonSchema


class MasterComponent(SkynetApiBase):
    __schema__ = MasterComponentSchema

    exp = relationship('expId', 'Exposure', uselist=False)
    master = relationship('masterId', 'MasterCalibration', uselist=False)


class MasterCalibration(SkynetApiBase):
    __schema__ = MasterCalibrationSchema
    __get_method__ = staticmethod(master_cal_methods.get)
    __query_method__ = staticmethod(master_cal_methods.query)

    components = relationship('components', 'MasterComponent')
    filter = relationship('filterId', 'StandardFilter')
    telescope = relationship('teleId', 'Telescope')

    @property
    def centerTime(self):
        return (
            self.date + timedelta(seconds=self.expLength/2)) if self.expLength \
            else self.date


class Observation(SkynetApiBase):
    __schema__ = ObservationSchema
    __get_method__ = staticmethod(obs_methods.get)
    __query_method__ = staticmethod(obs_methods.query)
    __update_method__ = staticmethod(obs_methods.update)

    user = relationship('userId', 'User', uselist=False)
    group = relationship('groupId', 'Group', uselist=False)
    collab = relationship('collabId', 'Collaboration', uselist=False)
    teleOwner = relationship('teleOwnerId', 'TeleOwner', uselist=False)
    userPreference = relationship('id', 'ObsUserPreference', uselist=False)
    currentTelescope = relationship('currentTeleId', 'Telescope', uselist=False)
    timeAccount = relationship('timeAccountId', 'TimeAccount', uselist=False)
    parentGroupTimeAccount = relationship(
        'parentGroupTimeAccountId', 'TimeAccount', uselist=False)
    parentCollabTimeAccount = relationship(
        'parentCollabTimeAccountId', 'TimeAccount', uselist=False)
    orbitalElements = relationship(
        'orbitalElementsId', 'OrbitalElements', uselist=False)
    trigger = relationship(
        'trigger', 'CampaignTrigger', uselist=False)

    obsNotes = relationship('obsNotes', 'ObsNote')
    exps = relationship('exps', 'Exposure')
    telescopes = relationship('telescopes', 'Telescope')

    def submit(self):
        """
        Submit an observation with parameters set by the current instance fields

        :return: None; the observation is updated in place
        :rtype: None
        """
        self.refresh(obs_methods.add(api_key=self._api_key, **self.fields))

    def activate(self):
        obs_methods.update(obs=self.id, state='active', api_key=self._api_key)
        self.refresh()

    def cancel(self):
        obs_methods.update(obs=self.id, state='canceled', api_key=self._api_key)
        self.refresh()


class ObsNote(SkynetApiBase):
    __schema__ = ObsNoteSchema


class ObsStats(SkynetApiBase):
    __schema__ = ObsStatsSchema


class ObsUserPreference(SkynetApiBase):
    __schema__ = ObsUserPreferenceSchema


class OrbitalElements(SkynetApiBase):
    __schema__ = OrbitalElementsSchema


class PendingCollabMembership(SkynetApiBase):
    __schema__ = PendingCollabMembershipSchema


class PendingGroupMembership(SkynetApiBase):
    __schema__ = PendingGroupMembershipSchema


class PriorityAccessSchedule(SkynetApiBase):
    __schema__ = PriorityAccessScheduleSchema


class RadioCartographerJob(DataProcessingJob):
    __schema__ = RadioCartographerJobSchema

    radioObs = relationship('radioObsId', 'RadioObs', uselist=False)


class RadioFilter(SkynetApiBase):
    __schema__ = RadioFilterSchema


class RadioObs(SkynetApiBase):
    __schema__ = RadioObsSchema
    __get_method__ = staticmethod(radio_obs_methods.get)
    __query_method__ = staticmethod(radio_obs_methods.query)
    __update_method__ = staticmethod(radio_obs_methods.update)
    __polymorphic_on__ = 'obsType'
    __polymorphic_table__ = {
        'track': 'RadioObsTrack',
        'map': 'RadioObsMap',
        'daisy': 'RadioObsDaisy',
        'onoff': 'RadioObsOnOff',
    }

    antennas = relationship('antennas', 'Antenna')
    currentAntenna = relationship('currentAntennaId', 'Antenna', uselist=False)
    radioFilter = relationship('radioFilterId', 'RadioFilter', uselist=False)
    receiver = relationship('receiverId', 'Receiver', uselist=False)
    user = relationship('userId', 'User', uselist=False)
    group = relationship('groupId', 'Group', uselist=False)
    collab = relationship('collabId', 'Collaboration', uselist=False)
    timeAccount = relationship('timeAccountId', 'TimeAccount', uselist=False)
    teleOwner = relationship('teleOwnerId', 'TeleOwner', uselist=False)
    leftJob = relationship('leftJob', 'RadioCartographerJob', uselist=False)
    rightJob = relationship('rightJob', 'RadioCartographerJob', uselist=False)
    compositeJob = relationship(
        'compositeJob', 'RadioCartographerJob', uselist=False)
    radioCartographerJobs = relationship(
        'radioCartographerJobs', 'RadioCartographerJob')

    def submit(self):
        """
        Submit a radio observation with parameters set by the current instance
        fields

        :return: None; the observation is updated in place
        :rtype: None
        """
        self.refresh(radio_obs_methods.add(
            api_key=self._api_key, **self.fields))

    def activate(self):
        """
        Activate the radio observation

        :return: None; the observation is updated in place
        :rtype: None
        """
        radio_obs_methods.update(
            obs=self.id, state='active', api_key=self._api_key)
        self.refresh()

    def cancel(self):
        """
        Cancel the radio observation

        :return: None; the observation is updated in place
        :rtype: None
        """
        radio_obs_methods.update(
            obs=self.id, state='canceled', api_key=self._api_key)
        self.refresh()

    def getBeamWidth(self, aperture_diameter=None):
        """
        Return the effective beam width

        :param float aperture_diameter: optional antenna aperture diameter in
            meters; if omitted, use the current antenna's diameter

        :return: beam width in degrees
        :rtype: float
        """
        if aperture_diameter is None:
            self.refresh()
            if self.currentAntenna is None:
                aperture_diameter = self.antennas[0].diameter
            else:
                aperture_diameter = self.currentAntenna.diameter

        if self.recordType in ('lowres', 'pulsar'):
            center_freq = self.centerFrequency
        elif self.recordType == 'highres':
            center_freq = max(self.centerFrequency, self.secondaryFrequency)
        else:
            raise Exception('Invalid receiver mode: {}'.format(self.recordType))

        center_wavelength_meters = 299.792458/center_freq

        return math.degrees(1.22*center_wavelength_meters/aperture_diameter)


class RadioObsTrack(RadioObs):
    __schema__ = RadioObsTrackSchema


class RadioObsMap(RadioObs):
    __schema__ = RadioObsMapSchema


class RadioObsDaisy(RadioObs):
    __schema__ = RadioObsDaisySchema


class RadioObsOnOff(RadioObs):
    __schema__ = RadioObsOnOffSchema


class ReceiverHistory(SkynetApiBase):
    __schema__ = ReceiverHistorySchema


class Receiver(SkynetApiBase):
    __schema__ = ReceiverSchema


class Role(SkynetApiBase):
    __schema__ = RoleSchema


class Site(SkynetApiBase):
    __schema__ = SiteSchema


class StandardFilter(Filter):
    __schema__ = StandardFilterSchema


class TeleOwnerGuestShare(SkynetApiBase):
    __schema__ = TeleOwnerGuestShareSchema


class TeleOwner(SkynetApiBase):
    __schema__ = TeleOwnerSchema


class TeleOwnerShare(SkynetApiBase):
    __schema__ = TeleOwnerShareSchema


class TelescopeBase(SkynetApiBase):
    __schema__ = TelescopeBaseSchema
    __get_method__ = staticmethod(tele_methods.get)
    __query_method__ = staticmethod(tele_methods.query)
    __polymorphic_on__ = 'teleType'
    __polymorphic_table__ = {'optical': 'Telescope', 'radio': 'Antenna'}

    site = relationship('siteId', 'Site', uselist=False)

    localHorizonData = relationship('localHorizonData', 'LocalHorizon')
    ownerShares = relationship('ownerShares', 'TeleOwnerShare')
    guestShares = relationship('guestShares', 'TeleOwnerGuestShare')
    status = relationship('id', 'TeleStatus', uselist=False)
    admins = relationship('admins', 'User')


class Telescope(TelescopeBase):
    __schema__ = TelescopeSchema

    dmsStatus = relationship('id', 'DmsStatus', uselist=False)
    filters = relationship('filters', 'StandardFilter')
    grbFilters = relationship('grbFilters', 'StandardFilter')
    darkCalibrationConfigs = relationship(
        'darkCalibrationConfigs', 'DarkCalibrationConfig')
    flatCalibrationConfigs = relationship(
        'flatCalibrationConfigs', 'FlatCalibrationConfig')
    # TODO: def cameraHistory(self):


class Antenna(TelescopeBase):
    __schema__ = AntennaSchema

    receiverHistory = relationship('receiverHistory', 'ReceiverHistory')
    receiver = relationship('receiverId', 'Receiver', uselist=False)


class TeleSiteRole(SkynetApiBase):
    __schema__ = TeleSiteRoleSchema


class TeleStatus2(SkynetApiBase):
    __schema__ = TeleStatus2Schema


class TimeAccount(SkynetApiBase):
    __schema__ = TimeAccountSchema

    @property
    def parentEntity(self):
        if self.parentGroup is not None:
            return self.parentGroup
        if self.parentCollab is not None:
            return self.parentCollab
        if self.parentTeleOwner is not None:
            return self.parentTeleOwner

        return None

    @property
    def parentTimeAccount(self):
        if self.parentGroupTimeAccount is not None:
            return self.parentGroupTimeAccount
        if self.parentCollabTimeAccount is not None:
            return self.parentCollabTimeAccount
        if self.parentTeleOwnerTimeAccount is not None:
            return self.parentTeleOwnerTimeAccount

        return None

    @property
    def ownerEntity(self):
        if self.user is not None:
            return self.user
        if self.group is not None:
            return self.group
        if self.collab is not None:
            return self.collab
        if self.teleOwner is not None:
            return self.teleOwner

        return None

    @property
    def ownerName(self):
        owner_entity = self.ownerEntity
        if owner_entity is None:
            return None
        if owner_entity is self.user:
            return owner_entity.username
        return owner_entity.name


class TimeAccountSettings(SkynetApiBase):
    __schema__ = TimeAccountSettingsSchema


class TimeAccountTxn(SkynetApiBase):
    __schema__ = TimeAccountTxnSchema


class UserSelfRegistration(SkynetApiBase):
    __schema__ = UserSelfRegistrationSchema

    user = relationship('user', 'User', uselist=False)


class User(SkynetApiBase):
    __schema__ = UserSchema
    __get_method__ = staticmethod(user_methods.get)
    __query_method__ = staticmethod(user_methods.query)
    __update_method__ = staticmethod(user_methods.update)

    selfRegistration = relationship(
        'selfRegistration', 'UserSelfRegistration', uselist=False)
    timeAccounts = relationship('timeAccounts', 'TimeAccount')
    groups = relationship('groups', 'Group')
    owningGroup = relationship('owningGroupId', 'Group', uselist=False)
    roles = relationship('roles', 'Role')
    userCollabRoles = relationship('userCollabRoles', 'Role')
    userGroupRoles = relationship('userGroupRoles', 'Role')
    userSiteRoles = relationship('userSiteRoles', 'Role')
    userTeleRoles = relationship('userTeleRoles', 'Role')
    userTeleOwnerRoles = relationship('userTeleOwnerRoles', 'Role')
    pendingGroupMemberships = relationship(
        'pendingGroupMemberships', 'PendingGroupMembership')

    # TODO: def skynetEventNotificationRequests(self):

    # TODO: def obsEventNotificationRequests(self):

    # TODO: def timeAccountEventNotificationRequests(self):

    # TODO: def messageEventNotificationRequests(self):

    # TODO: def teleEventNotificationRequests(self):

    # TODO: def siteEventNotificationRequests(self):

    @property
    def entityId(self):
        return 'u{:d}'.format(self.id)

    @property
    def entityType(self):
        return 'User'

    @property
    def name(self):
        return '{} {}'.format(self.firstName, self.lastName) \
            if self.firstName and self.lastName \
            else self.firstName if self.firstName \
            else self.lastName if self.lastName else ''

    @property
    def isMinor(self):
        return self.age is None or self.age < 18

    @property
    def age(self):
        if self.birthdate is None:
            return None
        bd = datetime.combine(self.birthdate, datetime.min.time())
        return (datetime.now() - bd).total_seconds()/(365.2425*86400)


class WcsSolution(SkynetApiBase):
    __schema__ = WcsSolutionSchema
