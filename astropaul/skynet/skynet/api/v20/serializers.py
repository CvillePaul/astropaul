"""
SkyNet API: JSON serializers for skynet.db.model classes
"""

from __future__ import absolute_import, division, print_function

import datetime
import json
from marshmallow import Schema, ValidationError, fields


# Treat empty string as None when assigning to non-string fields
def _serialize(self, value, attr, obj, **kwargs):
    """
    Serializer for non-string scalar fields that treats empty string as None on
    assignment to an API ORM object attribute

    :param marshmallow.fields.Field self: field object
    :param value: value to serialize
    :param str attr: schema attribute name
    :param obj: schema object

    :return: serialized value
    """
    if value is None or value == '':
        return None
    # noinspection PyProtectedMember
    return super(self.__class__, self)._serialize(value, attr, obj, **kwargs)


class Boolean(fields.Boolean):
    _serialize = _serialize


class Date(fields.Date):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for UTC date fields

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value is None or value == '':
            return None
        if isinstance(value, datetime.datetime):
            value = value.date()
        if not isinstance(value, datetime.date):
            try:
                value = self._deserialize(str(value), attr, obj, **kwargs)
            except (AttributeError, ValueError):
                self.fail('format', input=value)
        return value.isoformat()

    def _deserialize(self, value, attr, data, **kwargs):
        """
        Deserializer for date fields that accepts empty string, None, and
        datetime/date instance

        :param value: value to serialize; for datetime values, the date part is
            extracted
        :param str attr: schema attribute name
        :param data: unused

        :return: deserialized value
        """
        if value is None or value == '' and getattr(self, 'allow_none', False):
            return None
        try:
            return datetime.datetime.strptime(
                value.split('T', 1)[0], '%Y-%m-%d').date()
        except ValueError:
            raise self.fail('invalid')


class DateTime(fields.DateTime):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for UTC date/time fields

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value is None or value == '':
            return None
        if not isinstance(value, datetime.datetime):
            try:
                value = self._deserialize(str(value), attr, obj, **kwargs)
            except (AttributeError, ValueError):
                self.fail('format', input=value)
        return value.isoformat()

    def _deserialize(self, value, attr, data, **kwargs):
        """
        Deserializer for datetime fields that accepts empty string, None, and
        datetime/date instance

        :param value: value to serialize; for date-only values, midnight is
            assumed
        :param str attr: schema attribute name
        :param data: unused

        :return: deserialized value
        """
        if value is None or value == '' and getattr(self, 'allow_none', False):
            return None
        try:
            return datetime.datetime.strptime(value, '%Y-%m-%dT%H:%M:%S.%f')
        except ValueError:
            try:
                return datetime.datetime.strptime(value, '%Y-%m-%dT%H:%M:%S')
            except ValueError:
                raise self.fail('invalid')


class Email(fields.Email):
    def _serialize(self, value, *args, **kwargs):
        """
        Serializer for email fields that appends "INVALID" to non-conforming
        emails

        :param value: value to serialize
        :param str attr: schema attribute name
        :param marshmallow.Schema obj: schema object

        :return: serialized value
        """
        if value == '' or value is None:
            return None
        try:
            return super(self.__class__, self)._serialize(
                value, *args, **kwargs)
        except ValidationError:
            return '<INVALID> {}'.format(value)


class Float(fields.Float):
    _serialize = _serialize


class FloatSex(fields.Float):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for floating-point fields that store angles, e.g. RA and Dec;
        accepts "+DD:MM:SS.S" on assignment

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value == '' or value is None:
            return None
        if isinstance(value, str) or isinstance(value, type(u'')):
            try:
                d, m, s = value.split(':')
                value = (abs(int(d)) + int(m)/60 + float(s)/3600) * \
                    (1 - 2*(d.strip()[:1] == '-'))
            except ValueError:
                try:
                    d, m = value.split(':')
                    value = (abs(int(d)) + float(m)/60) * \
                        (1 - 2*(d.strip()[:1] == '-'))
                except ValueError:
                    pass
        return super(self.__class__, self)._serialize(
            value, attr, obj, **kwargs)


class Integer(fields.Integer):
    _serialize = _serialize


class List(fields.List):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for list fields that treats empty string as None on
        assignment to an API ORM object attribute and allows assignments from
        comma-separated strings of values; this is intended for an easier use on
        the client side in situations like
            obs = Observation(..., telescopes='Prompt1,Prompt2')
        instead of
            obs = Observation(...,
                telescopes=[{'name': 'Prompt1'}, {'name': 'Prompt2'}])

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value == '' or value is None:
            return None

        try:
            # noinspection PyUnresolvedReferences
            inner = self.container  # marshmallow 2
        except AttributeError:
            # noinspection PyUnresolvedReferences
            inner = self.inner  # marshmallow 3
        if isinstance(inner, fields.Nested) and (
                isinstance(value, str) or isinstance(value, type(u''))):
            if value.startswith('[') and value.endswith(']'):
                # A JSON list
                value = json.loads(value)
            else:
                # A comma-separated list of values
                value = [s.strip() for s in value.split(',')]

        return super(List, self)._serialize(value, attr, obj, **kwargs)


class Nested(fields.Nested):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for nested fields that treats empty string as None on
        assignment to an API ORM object attribute and allows for assignments
        like field = 1, field = '1' (both are the same as field = {'id': 1}),
        and field = 'name' (same as field = {'name': 'bar'}), treating integers
        as nested object IDs and non-integers as names and also accepting JSON
        strings

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value == '' or value is None:
            return None

        # Convert scalars to nested schema representations with either ID or
        # name
        if isinstance(value, str) or isinstance(value, type(u'')):
            try:
                decoded_value = json.loads(value)
                if not isinstance(decoded_value, dict):
                    raise ValueError()

                # Value is a JSON structure string, but does it match the
                # current schema? Exclude the (possibly extremely rare)
                # cases when the nested field name looks like a JSON
                # structure by checking that at least one field name in the
                # value being assigned is the name of one of the nested schema
                # fields.
                # noinspection PyProtectedMember
                if not any(key in self.nested._declared_fields
                           for key in decoded_value):
                    raise ValueError()
            except ValueError:
                # Treat other strings as IDs or names
                id_field = getattr(self.schema, '__id_field__', 'id')
                name_field = getattr(self.schema, '__name_field__', 'name')
                if not id_field and not name_field:
                    raise ValueError(
                        'Defining "{}" objects by ID/name is not '
                        'supported'.format(self.schema.__class__.__name__))
                if id_field:
                    try:
                        value = {id_field: int(value)}
                    except ValueError:
                        if name_field:
                            value = {name_field: value}
                        else:
                            raise ValueError(
                                'Invalid "{}" object ID'.format(
                                    self.schema.__class__.__name__))
                else:
                    value = {name_field: value}
            else:
                # Treat as JSON string containing the current nested field
                value = decoded_value

        elif isinstance(value, int):
            id_field = getattr(self.schema, '__id_field__', 'id')
            if not id_field:
                raise ValueError(
                    'Defining "{}" objects by ID is not supported'.format(
                        self.schema.__class__.__name__))
            value = {id_field: value}

        return super(Nested, self)._serialize(value, attr, obj, **kwargs)


class Time(fields.Time):
    def _serialize(self, value, attr, obj, **kwargs):
        """
        Serializer for UTC time fields

        :param value: value to serialize
        :param str attr: schema attribute name
        :param obj: schema object

        :return: serialized value
        """
        if value is None or value == '':
            return None
        if isinstance(value, datetime.datetime):
            value = value.time()
        if not isinstance(value, datetime.time):
            try:
                value = self._deserialize(str(value), attr, obj, **kwargs)
            except (AttributeError, ValueError):
                self.fail('format', input=value)
        return value.isoformat()

    def _deserialize(self, value, attr, data, **kwargs):
        """
        Deserializer for time fields that accepts empty string, None, and
        datetime/time instance

        :param value: value to serialize; for datetime values, the time part is
            extracted
        :param str attr: schema attribute name
        :param data: unused

        :return: deserialized value
        """
        if value is None or value == '' and getattr(self, 'allow_none', False):
            return None
        value = value.split('T', 1)[-1]
        try:
            return datetime.datetime.strptime(value, '%H:%M:%S.%f').time()
        except ValueError:
            try:
                return datetime.datetime.strptime(value, '%H:%M:%S').time()
            except ValueError:
                raise self.fail('invalid')


String = fields.String


def is_scalar_field(field):
    """
    Is the given field scalar (i.e. directly stored in a database table, e.g.
    Integeer and String) or composite (i.e. corresponds to a relationship in
    the corresponding db.model class)?

    :param marshmallow.fields.Field field: schema field

    :return: True for scalar fields
    :rtype: bool
    """
    return not any(
        isinstance(field, t)
        for t in [
            fields.Nested, fields.List, fields.Dict, fields.Method,
            fields.Function])


# By default, serializers will marshal the object attributes that have the same
# name as the fields. Included attribute names for clarity, but note that they
# are not *necessary*. This will also allow us to change the attribute names in
# the future easily.

class SkynetSchema(Schema):
    def __init__(self, *args, **kwargs):
        # Set Schema.partial = True and Schema.unknown = 'include' by
        # default -- required for API ORM to function properly
        kwargs = dict(kwargs)
        if 'partial' not in kwargs:
            kwargs['partial'] = True
        if 'unknown' not in kwargs:
            kwargs['unknown'] = 'include'
        super().__init__(*args, **kwargs)


class CameraSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')
    pixelWidth = Integer(attribute='pixelWidth')
    pixelHeight = Integer(attribute='pixelHeight')
    pixelSize = Float(attribute='pixelSize', allow_none=True)
    blankLevel = Float(attribute='blankLevel', dump_only=True)
    gain = Float(attribute='gain', dump_only=True)
    readNoise = Float(attribute='readNoise', dump_only=True)


class CampaignLimitedSchema(SkynetSchema):
    __name_field__ = None  # campaigns are identified by ID, name is not unique

    id = Integer(attribute='id')
    name = String(attribute='name')


class CampaignTriggerLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    campaignId = Integer(attribute='campaignId')


class CollaborationLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class DmsStatusSchema(SkynetSchema):
    __id_field__ = 'teleId'
    __name_field__ = None

    teleId = Integer(attribute='teleId')
    expId = Integer(attribute='expId', allow_none=True)
    fileSize = Integer(attribute='fileSize', allow_none=True)
    bytesTransferred = Integer(
        attribute='bytesTransferred', allow_none=True)
    kilobytesPerSecond = Float(
        attribute='kilobytesPerSecond', allow_none=True)
    filename = String(attribute='filename', allow_none=True)
    lastUpdate = DateTime(attribute='lastUpdate')


class ExposureLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    obsId = Integer(attribute='obsId')


class FilterLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class FilterSchema(FilterLimitedSchema):
    order = Integer(attribute='order', allow_none=True)
    type = String(attribute='type')
    throughput = Float(attribute='throughput')
    flatScaler = Float(attribute='flatScaler')


class GroupLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class LocalHorizonSchema(SkynetSchema):
    __id_field__ = None
    __name_field__ = None

    azDegs = FloatSex(attribute='azDegs')
    elDegs = FloatSex(attribute='elDegs')


class ObservationLimitedSchema(SkynetSchema):
    __name_field__ = None  # obs are identified by ID, name is not unique

    id = Integer(attribute='id')
    name = String(attribute='name')


class ObsNoteSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    obsId = Integer(attribute='obsId')
    note = String(attribute='note')
    timeCreated = DateTime(attribute='timeCreated')

    obs = Nested(ObservationLimitedSchema, attribute='obs', allow_none=True)


class ObsStatsSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    radioObsCompleted = Integer(attribute='radioObsCompleted')
    radioObsArchived = Integer(attribute='radioObsArchived')
    opticalExpCompleted = Integer(attribute='opticalExpCompleted')
    opticalExpArchived = Integer(attribute='opticalExpArchived')


class ObsUserPreferenceLimitedSchema(SkynetSchema):
    __id_field__ = 'obsId'
    __name_field__ = None

    obsId = Integer(attribute='obsId')
    userId = Integer(attribute='userId')
    reduce = Boolean(attribute='reduce')
    invert = Boolean(attribute='invert')
    scalePreset = String(attribute='scalePreset')


class ObsUserPreferenceSchema(ObsUserPreferenceLimitedSchema):
    obs = Nested(ObservationLimitedSchema, attribute='obs', allow_none=True)


class OrbitalElementsSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    epochJd = Float(attribute='epochJd')
    m = Float(attribute='m')
    peri = Float(attribute='peri')
    node = Float(attribute='node')
    incl = Float(attribute='incl')
    e = Float(attribute='e')
    a = Float(attribute='a')


class PriorityAccessScheduleLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')


class RadioObsLimitedSchema(SkynetSchema):
    __name_field__ = None  # radio obs are identified by ID, name is not unique

    id = Integer(attribute='id')
    name = String(attribute='name')
    obsType = String(attribute='obsType')


class RadioFilterLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class RadioFilterSchema(RadioFilterLimitedSchema):
    description = String(attribute='description')
    startFreq = Float(attribute='startFreq', allow_none=True)
    stopFreq = Float(attribute='stopFreq', allow_none=True)


class ReceiverLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class RoleSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')
    description = String(attribute='description')


class SiteLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')
    latDegs = FloatSex(attribute='latDegs')
    lngDegs = FloatSex(attribute='lngDegs')


class TeleOwnerGuestShareLimitedSchema(SkynetSchema):
    __id_field__ = 'teleId'
    __name_field__ = None

    teleId = Integer(attribute='teleId')
    teleOwnerId = Integer(attribute='teleOwnerId')
    totalShare = Float(attribute='totalShare')
    contesting = Integer(attribute='contesting')
    totalUsage = Float(attribute='totalUsage')
    totalContestedUsage = Float(attribute='totalContestedUsage')
    totalTimeWaiting = Float(attribute='totalTimeWaiting')
    pendingTotalUsage = Float(attribute='pendingTotalUsage')
    pendingTotalContestedUsage = Float(attribute='pendingTotalContestedUsage')
    pendingTotalTimeWaiting = Float(attribute='pendingTotalTimeWaiting')
    lastUpdate = DateTime(attribute='lastUpdate', allow_none=True)


class TeleOwnerLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')


class TeleOwnerShareLimitedSchema(SkynetSchema):
    __id_field__ = 'teleId'
    __name_field__ = None

    teleId = Integer(attribute='teleId')
    teleOwnerId = Integer(attribute='teleOwnerId')
    totalShare = Float(attribute='totalShare')
    ownerShare = Float(attribute='ownerShare')
    guestShare = Float(attribute='guestShare')
    publicShare = Float(attribute='publicShare')
    ownerContesting = Integer(attribute='ownerContesting')
    guestContesting = Integer(attribute='guestContesting')
    publicContesting = Integer(attribute='publicContesting')
    hasTooPrivilege = Boolean(attribute='hasTooPrivilege')
    totalUsage = Float(attribute='totalUsage')
    totalContestedUsage = Float(attribute='totalContestedUsage')
    totalTimeWaiting = Float(attribute='totalTimeWaiting')
    ownerUsage = Float(attribute='ownerUsage')
    ownerContestedUsage = Float(attribute='ownerContestedUsage')
    ownerTimeWaiting = Float(attribute='ownerTimeWaiting')
    guestUsage = Float(attribute='guestUsage')
    guestContestedUsage = Float(attribute='guestContestedUsage')
    guestTimeWaiting = Float(attribute='guestTimeWaiting')
    publicUsage = Float(attribute='publicUsage')
    publicContestedUsage = Float(attribute='publicContestedUsage')
    publicTimeWaiting = Float(attribute='publicTimeWaiting')
    pendingTotalUsage = Float(attribute='pendingTotalUsage')
    pendingTotalContestedUsage = Float(attribute='pendingTotalContestedUsage')
    pendingTotalTimeWaiting = Float(attribute='pendingTotalTimeWaiting')
    pendingOwnerUsage = Float(attribute='pendingOwnerUsage')
    pendingOwnerContestedUsage = Float(attribute='pendingOwnerContestedUsage')
    pendingOwnerTimeWaiting = Float(attribute='pendingOwnerTimeWaiting')
    pendingGuestUsage = Float(attribute='pendingGuestUsage')
    pendingGuestContestedUsage = Float(attribute='pendingGuestContestedUsage')
    pendingGuestTimeWaiting = Float(attribute='pendingGuestTimeWaiting')
    pendingPublicUsage = Float(attribute='pendingPublicUsage')
    pendingPublicContestedUsage = Float(attribute='pendingPublicContestedUsage')
    pendingPublicTimeWaiting = Float(attribute='pendingPublicTimeWaiting')
    lastUpdate = DateTime(attribute='lastUpdate', allow_none=True)


class TeleSiteRoleSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')
    description = String(attribute='description')
    siteOnly = Boolean(attribute='siteOnly')


class TeleStatus2LimitedSchema(SkynetSchema):
    __id_field__ = 'teleId'
    __name_field__ = None

    teleId = Integer(attribute='teleId')
    connected = Boolean(attribute='connected')
    cameraStatus = String(attribute='cameraStatus', allow_none=True)
    domeStatus = String(attribute='domeStatus', allow_none=True)
    mountStatus = String(attribute='mountStatus', allow_none=True)
    weatherStatus = String(attribute='weatherStatus', allow_none=True)
    obsId = Integer(attribute='obsId', allow_none=True)
    radioObsId = Integer(attribute='radioObsId', allow_none=True)
    lastUpdate = DateTime(attribute='lastUpdate', allow_none=True)
    raHours = FloatSex(attribute='raHours', allow_none=True)
    decDegs = FloatSex(attribute='decDegs', allow_none=True)
    azDegs = FloatSex(attribute='azDegs', allow_none=True)
    elDegs = FloatSex(attribute='elDegs', allow_none=True)
    expProgress = Float(attribute='expProgress', allow_none=True)
    expLength = Float(attribute='expLength', allow_none=True)
    expId = Integer(attribute='expId', allow_none=True)
    focusStatus = String(attribute='focusStatus', allow_none=True)
    operatingMode = String(attribute='operatingMode', allow_none=True)


class TelescopeBaseLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')
    name = String(attribute='name')
    site = Nested(SiteLimitedSchema, attribute='site', allow_none=True)
    siteName = String(attribute='siteName')


class TimeAccountLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')


class UserLimitedSchema(SkynetSchema):
    __name_field__ = 'username'

    id = Integer(attribute='id')
    username = String(attribute='username', allow_none=True)


class WcsSolutionLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')


class PendingGroupMembershipLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')


class PendingGroupMembershipSchema(PendingGroupMembershipLimitedSchema):
    groupId = Integer(attribute='groupId')
    userId = Integer(attribute='userId')
    createdOn = DateTime(attribute='createdOn')
    state = String(attribute='state')
    isInvitation = Boolean(attribute='isInvitation')

    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    user = Nested(UserLimitedSchema, attribute='user', allow_none=True)


class PendingCollabMembershipLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')


class PendingCollabMembershipSchema(PendingCollabMembershipLimitedSchema):
    groupId = Integer(attribute='groupId')
    collabId = Integer(attribute='collabId')
    createdOn = DateTime(attribute='createdOn')
    state = String(attribute='state')
    isInvitation = Boolean(attribute='isInvitation')

    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    collab = Nested(
        CollaborationLimitedSchema, attribute='collab', allow_none=True)


class UserSelfRegistrationSchema(SkynetSchema):
    id = Integer(attribute='id')
    keyCreatedOn = DateTime(attribute='keyCreatedOn')
    registrationKey = String(attribute='registrationKey')
    activated = Boolean(attribute='activated')
    user = Nested(UserLimitedSchema, attribute='user')


class UserAdminSchema(UserLimitedSchema):
    owningGroupId = Integer(attribute='owningGroupId', allow_none=True)
    priority = Integer(attribute='priority')
    level = Integer(attribute='level')
    canModify = Boolean(attribute='canModify')
    groupDeactivated = Boolean(attribute='groupDeactivated')
    deleted = Boolean(attribute='deleted', allow_none=True)


class UserRestrictedSchema(UserAdminSchema):
    alertsOn = Integer(attribute='alertsOn')
    reduceImagesByDefault = Boolean(attribute='reduceImagesByDefault')
    obsViewMode = String(attribute='obsViewMode')
    expViewMode = String(attribute='expViewMode')
    prefer16BitImages = Boolean(attribute='prefer16BitImages')
    hideClosedOwnerTimeAccounts = Boolean(
        attribute='hideClosedOwnerTimeAccounts')
    hideClosedCollabTimeAccounts = Boolean(
        attribute='hideClosedCollabTimeAccounts')
    hideClosedGroupTimeAccounts = Boolean(
        attribute='hideClosedGroupTimeAccounts')
    apiExpCheckMode = String(attribute='apiExpCheckMode')


class TimeAccountSettingsSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    maxPriority = Integer(attribute='maxPriority')
    teleOwnerId = Integer(attribute='teleOwnerId', allow_none=True)
    accountType = String(attribute='accountType')
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    telescopes = List(
        Nested(TelescopeBaseLimitedSchema), attribute='telescopes',
        dump_only=True)


class TimeAccountSchema(TimeAccountLimitedSchema):
    teleOwnerId = Integer(attribute='teleOwnerId', allow_none=True)
    collabId = Integer(attribute='collabId', allow_none=True)
    groupId = Integer(attribute='groupId', allow_none=True)
    userId = Integer(attribute='userId', allow_none=True)
    parentTeleOwnerTimeAccountId = Integer(
        attribute='parentTeleOwnerTimeAccountId', allow_none=True)
    parentTeleOwnerId = Integer(attribute='parentTeleOwnerId', allow_none=True)
    parentCollabTimeAccountId = Integer(
        attribute='parentCollabTimeAccountId', allow_none=True)
    parentCollabId = Integer(attribute='parentCollabId', allow_none=True)
    parentGroupTimeAccountId = Integer(
        attribute='parentGroupTimeAccountId', allow_none=True)
    parentGroupId = Integer(attribute='parentGroupId', allow_none=True)
    timeAccountSettingsId = Integer(attribute='timeAccountSettingsId')
    description = String(attribute='description', allow_none=True)
    active = Boolean(attribute='active')

    parentTeleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='parentTeleOwner', allow_none=True)
    parentCollab = Nested(
        CollaborationLimitedSchema, attribute='parentCollab', allow_none=True)
    parentGroup = Nested(
        GroupLimitedSchema, attribute='parentGroup', allow_none=True)

    # txns = List(Nested(
    #     TimeAccountTxnSchema), attribute='txns')
    settings = Nested(
        TimeAccountSettingsSchema, attribute='settings', allow_none=True)
    balance = Integer(attribute='balance', dump_only=True)
    creditsReceived = Integer(attribute='creditsReceived', dump_only=True)
    creditsUsed = Integer(attribute='creditsUsed', dump_only=True)

    user = Nested(UserLimitedSchema, attribute='user', allow_none=True)
    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    collab = Nested(
        CollaborationLimitedSchema, attribute='collab', allow_none=True)
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)


class UserSchema(UserRestrictedSchema):
    email = Email(attribute='email', allow_none=True)
    firstName = String(attribute='firstName', allow_none=True)
    lastName = String(attribute='lastName', allow_none=True)
    phoneEmail = String(attribute='phoneEmail', allow_none=True)
    birthdate = Date(attribute='birthdate', allow_none=True)
    consumerKey = String(attribute='consumerKey', allow_none=True)
    consumerSecret = String(attribute='consumerSecret', allow_none=True)
    accessToken = String(attribute='accessToken', allow_none=True)
    location = String(attribute='location', allow_none=True)
    userSelfRegistrationId = Integer(
        attribute='userSelfRegistrationId', allow_none=True)
    selfRegistration = Nested(
        UserSelfRegistrationSchema, attribute='selfRegistration',
        allow_none=True)

    timeAccounts = List(Nested(TimeAccountSchema), attribute='timeAccounts')
    # oauthClients = relationship(
    #     u'OauthClient', secondary='users_oauth_clients')

    groups = List(
        Nested(GroupLimitedSchema), attribute='groups')
    owningGroup = Nested(
        GroupLimitedSchema, attribute='owningGroup', allow_none=True)
    roles = List(Nested(RoleSchema), attribute='roles')
    userCollabRoles = List(Nested(RoleSchema), attribute='userCollabRoles')
    userGroupRoles = List(Nested(RoleSchema), attribute='userGroupRoles')
    userSiteRoles = List(Nested(RoleSchema), attribute='userSiteRoles')
    userTeleRoles = List(Nested(RoleSchema), attribute='userTeleRoles')
    userTeleOwnerRoles = List(
        Nested(RoleSchema), attribute='userTeleOwnerRoles')

    pendingGroupMemberships = List(
        Nested(PendingGroupMembershipLimitedSchema),
        attribute='pendingGroupMemberships')

    # TODO: skynetEventNotificationRequests = relationship(
    #     u'SkynetEventNotificationRequest', backref='user')
    # obsEventNotificationRequests = relationship(
    #     u'ObsEventNotificationRequest', backref='user')
    # timeAccountEventNotificationRequests = relationship(
    #     u'TimeAccountEventNotificationRequest', backref='user')
    # messageEventNotificationRequests = relationship(
    #     u'MessageEventNotificationRequest', backref='user')
    # teleEventNotificationRequests = relationship(
    #     u'TeleEventNotificationRequest', backref='user')
    # siteEventNotificationRequests = relationship(
    #     u'SiteEventNotificationRequest', backref='user')

    isSkynetAdmin = Boolean(attribute='isSkynetAdmin', dump_only=True)
    name = String(attribute='name', dump_only=True)


class WcsSolutionSchema(WcsSolutionLimitedSchema):
    expId = Integer(attribute='expId')
    foundSolution = Boolean(attribute='foundSolution')
    solutionType = Integer(attribute='solutionType')
    rejected = Boolean(attribute='rejected')
    crpix1 = Float(attribute='crpix1', allow_none=True)
    crpix2 = Float(attribute='crpix2', allow_none=True)
    crval1 = Float(attribute='crval1', allow_none=True)
    crval2 = Float(attribute='crval2', allow_none=True)
    cd11 = Float(attribute='cd11', allow_none=True)
    cd12 = Float(attribute='cd12', allow_none=True)
    cd21 = Float(attribute='cd21', allow_none=True)
    cd22 = Float(attribute='cd22', allow_none=True)
    width = Integer(attribute='width', allow_none=True)
    height = Integer(attribute='height', allow_none=True)
    rotation = Float(attribute='rotation', allow_none=True)
    pixelScale = Float(attribute='pixelScale', allow_none=True)
    mirrored = Boolean(attribute='mirrored', allow_none=True)
    dateSolved = DateTime(attribute='dateSolved', allow_none=True)
    pointingError = Float(attribute='pointingError', allow_none=True)
    deltaRa = Float(attribute='deltaRa', allow_none=True)
    deltaDec = Float(attribute='deltaDec', allow_none=True)
    logOdds = Float(attribute='logOdds', allow_none=True)
    nMatch = Integer(attribute='nMatch', allow_none=True)
    nConflict = Integer(attribute='nConflict', allow_none=True)
    nField = Integer(attribute='nField', allow_none=True)
    secondsToSolve = Float(attribute='secondsToSolve', allow_none=True)
    astrometryDotNetConfig = String(
        attribute='astrometryDotNetConfig', allow_none=True)
    waitTime = Float(attribute='waitTime', allow_none=True)
    dateProcessed = DateTime(attribute='dateProcessed', allow_none=True)
    numSources = Integer(attribute='numSources')
    fwhm = Float(attribute='fwhm', allow_none=True)
    posAngle = Float(attribute='posAngle', allow_none=True)
    ellipticity = Float(attribute='ellipticity', allow_none=True)
    zeroPoint = Float(attribute='zeroPoint', allow_none=True)
    skyBrightness = Float(attribute='skyBrightness', allow_none=True)
    limitingMag = Float(attribute='limitingMag', allow_none=True)
    limitingMagUnitExp = Float(attribute='limitingMagUnitExp', allow_none=True)
    efficiency = Float(attribute='efficiency', allow_none=True)


class ExposureSchema(ExposureLimitedSchema):
    expNum = Integer(attribute='expNum')
    expLength = Float(attribute='expLength')
    expLengthUsed = Float(attribute='expLengthUsed', allow_none=True)
    state = String(attribute='state')
    type = String(attribute='type')
    timeIn = DateTime(attribute='timeIn')
    startAfter = DateTime(attribute='startAfter', allow_none=True)
    endBefore = DateTime(attribute='endBefore', allow_none=True)
    timeTaken = DateTime(attribute='timeTaken', allow_none=True)
    timeSubmitted = DateTime(attribute='timeSubmitted', allow_none=True)
    timeArchived = DateTime(attribute='timeArchived', allow_none=True)
    isCompressed = Boolean(attribute='isCompressed')
    commandedRa = Float(attribute='commandedRa', allow_none=True)
    commandedDec = Float(attribute='commandedDec', allow_none=True)
    delay = Integer(attribute='delay', allow_none=True)
    timeTakenIsFromHdr = Boolean(attribute='timeTakenIsFromHdr')
    compression = String(attribute='compression', allow_none=True)
    creditsCharged = Integer(attribute='creditsCharged', allow_none=True)
    teleOwnerIdUsed = Integer(attribute='teleOwnerIdUsed', allow_none=True)
    timeAccountId = Integer(attribute='timeAccountId', allow_none=True)
    parentCollabTimeAccountId = Integer(
        attribute='parentCollabTimeAccountId', allow_none=True)
    parentGroupTimeAccountId = Integer(
        attribute='parentGroupTimeAccountId', allow_none=True)
    teleId = Integer(attribute='teleId', allow_none=True)
    teleIdRequested = Integer(attribute='teleIdRequested', allow_none=True)
    filterIdRequested = Integer(attribute='filterIdRequested', allow_none=True)
    filterIdUsed = Integer(attribute='filterIdUsed', allow_none=True)
    wcsId = Integer(attribute='wcsId', allow_none=True)
    wcsState = String(attribute='wcsState')
    statsState = String(attribute='statsState')
    numFileErrors = Integer(attribute='numFileErrors')
    targetExpId = Integer(attribute='targetExpId', allow_none=True)
    binningRequested = Integer(attribute='binningRequested', allow_none=True)
    binningUsed = Integer(attribute='binningUsed')

    linkedExps = List(Nested(ExposureLimitedSchema), attribute='linkedExps')

    filterRequested = Nested(
        FilterLimitedSchema, attribute='filterRequested', allow_none=True)
    filterUsed = Nested(
        FilterLimitedSchema, attribute='filterUsed', allow_none=True)
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)
    telescopeRequested = Nested(
        TelescopeBaseLimitedSchema, attribute='telescopeRequested',
        allow_none=True)
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    wcsSolution = Nested(
        WcsSolutionSchema, attribute='wcsSolution', allow_none=True)

    obs = Nested(ObservationLimitedSchema, attribute='obs', allow_none=True)

    centerTime = DateTime(attribute='centerTime', dump_only=True)
    filterName = String(attribute='filterName', dump_only=True)
    jpgDownloadUrl = String(attribute='jpgDownloadUrl', dump_only=True)
    fitsDownloadUrl = String(attribute='fitsDownloadUrl', dump_only=True)
    remoteFilename = String(attribute='remoteFilename', dump_only=True)
    camera = Nested(CameraSchema, attribute='camera', dump_only=True)
    blankLevel = Float(attribute='blankLevel', dump_only=True)
    gain = Float(attribute='gain', dump_only=True)
    readNoise = Float(attribute='readNoise', dump_only=True)
    imageId = String(attribute='imageId', dump_only=True)


class MasterComponentSchema(SkynetSchema):
    expId = Integer(attribute='expId')
    masterId = Integer(attribute='masterId')
    included = Boolean(attribute='included')
    rejectionReason = String(attribute='rejectionReason', allow_none=True)
    rejectionDescription = String(
        attribute='rejectionDescription', allow_none=True)
    biasMasterId = Integer(attribute='biasMasterId', allow_none=True)
    darkMasterId = Integer(attribute='darkMasterId', allow_none=True)

    exp = Nested(ExposureLimitedSchema)


class MasterCalibrationSchema(SkynetSchema):
    id = Integer(attribute='id')
    date = DateTime(attribute='date')
    startDate = DateTime('startDate')
    stopDate = DateTime('stopDate')
    combineMethod = String(attribute='combineMethod', default='mean')
    rejectMethod = String(attribute='rejectMethod', default='none')
    average = Float(attribute='average')
    stdev = Float(attribute='stdev')
    expLength = Float(attribute='expLength', allow_none=True)
    filterName = String(attribute='filterName', allow_none=True)
    imageType = String(attribute='imageType')
    scaleMethod = String(attribute='scaleMethod', default='none')
    rejected = Boolean(attribute='rejected', allow_none=True)
    rejectedBy = String(attribute='rejectedBy', allow_none=True)
    chauvenetRejectedPercentage = Float(
        attribute='chauvenetRejectedPercentage', allow_none=True)
    compression = String(attribute='compression', allow_none=True)
    onNewRaid = Boolean(attribute='onNewRaid', default=False)
    teleId = Integer(attribute='teleId')
    filterId = Integer(attribute='filterId', allow_none=True)
    binning = Integer(attribute='binning')

    components = List(Nested(MasterComponentSchema))
    filter = Nested(FilterLimitedSchema, allow_none=True)
    telescope = Nested(TelescopeBaseLimitedSchema)


class CampaignSchema(CampaignLimitedSchema):
    timeAccountId = Integer(attribute='timeAccountId')

    telescopes = List(
        Nested(TelescopeBaseLimitedSchema), attribute='telescopes')


class CampaignTriggerSchema(CampaignTriggerLimitedSchema):
    obsId = Integer(attribute='obsId')
    state = String(attribute='state')
    startTime = DateTime(attribute='startTime', allow_none=True)
    eventTime = DateTime(attribute='eventTime', allow_none=True)
    refTime = Float(attribute='refTime')
    refFilterId = Integer(attribute='refFilterId')
    refMagnitude = Float(attribute='refMagnitude')
    refSNR = Float(attribute='refSNR')
    tempIndex = Float(attribute='tempIndex')
    specIndex = Float(attribute='specIndex')
    expFunction = String(attribute='expFunction')
    maxCampaignLength = Float(attribute='maxCampaignLength', allow_none=True)

    obs = Nested(ObservationLimitedSchema, attribute='obs')
    campaign = Nested(CampaignLimitedSchema, attribute='campaign')
    refFilter = Nested(FilterLimitedSchema)
    refFilterRequested = Nested(
        FilterLimitedSchema, attribute='refFilterRequested', allow_none=True)


class ObservationSchema(ObservationLimitedSchema):
    raHours = FloatSex(attribute='raHours')
    decDegs = FloatSex(attribute='decDegs')
    grbId = Integer(attribute='grbId')
    state = String(attribute='state')
    minEl = Float(attribute='minEl')
    maxSun = Float(attribute='maxSun')
    mode = Integer(attribute='mode')
    timeIn = DateTime(attribute='timeIn')
    type = String(attribute='type')
    objectName = String(attribute='objectName', allow_none=True)
    objectType = String(attribute='objectType')
    objectDist = Float(attribute='objectDist', allow_none=True)
    cancelAfterUtc = DateTime(attribute='cancelAfterUtc', allow_none=True)
    timeSubmittedUtc = DateTime(attribute='timeSubmittedUtc', allow_none=True)
    totalTimeRemaining = Float(attribute='totalTimeRemaining')
    userId = Integer(attribute='userId', allow_none=True)
    nextExpStartAfterUtc = DateTime(
        attribute='nextExpStartAfterUtc', allow_none=True)
    collabId = Integer(attribute='collabId', allow_none=True)
    groupId = Integer(attribute='groupId', allow_none=True)
    timeAccountId = Integer(attribute='timeAccountId', allow_none=True)
    parentCollabTimeAccountId = Integer(
        attribute='parentCollabTimeAccountId', allow_none=True)
    parentGroupTimeAccountId = Integer(
        attribute='parentGroupTimeAccountId', allow_none=True)
    teleOwnerId = Integer(attribute='teleOwnerId', allow_none=True)
    accountType = String(attribute='accountType', allow_none=True)
    ditherEnabled = Boolean(attribute='ditherEnabled')
    ditherXSize = Integer(attribute='ditherXSize', allow_none=True)
    ditherYSize = Integer(attribute='ditherYSize', allow_none=True)
    ditherSpacingArcsecs = Float(
        attribute='ditherSpacingArcsecs', allow_none=True)
    targetTracking = String(attribute='targetTracking')
    triggerRepointEnabled = Boolean(attribute='triggerRepointEnabled')
    triggerRepointArcmins = Float(
        attribute='triggerRepointArcmins', allow_none=True)
    pointAheadEnabled = Boolean(attribute='pointAheadEnabled', allow_none=True)
    pointAheadSecs = Float(attribute='pointAheadSecs', allow_none=True)
    constantRaOffsetArcmins = Float(
        attribute='constantRaOffsetArcmins', allow_none=True)
    constantDecOffsetArcmins = Float(
        attribute='constantDecOffsetArcmins', allow_none=True)
    minMoonSepDegs = FloatSex(attribute='minMoonSepDegs', allow_none=True)
    fieldLockUtc = DateTime(attribute='fieldLockUtc', allow_none=True)
    rbiFractionAvgBkgLimit = Float(
        attribute='rbiFractionAvgBkgLimit', allow_none=True)
    tooJustification = String(attribute='tooJustification', allow_none=True)
    currentTeleId = Integer(attribute='currentTeleId', allow_none=True)
    priority = Integer(attribute='priority')
    efficiency = Float(attribute='efficiency', allow_none=True)
    orbitalElementsId = Integer(attribute='orbitalElementsId', allow_none=True)
    isToo = Boolean(attribute='isToo')

    user = Nested(UserLimitedSchema, attribute='user', allow_none=True)
    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    collab = Nested(
        CollaborationLimitedSchema, attribute='collab', allow_none=True)
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    userPreference = Nested(
        ObsUserPreferenceLimitedSchema, attribute='userPreference',
        allow_none=True)
    currentTelescope = Nested(
        TelescopeBaseLimitedSchema, attribute='currentTelescope',
        allow_none=True)
    orbitalElements = Nested(
        OrbitalElementsSchema, attribute='orbitalElements', allow_none=True)
    trigger = Nested(
        CampaignTriggerSchema, attribute='trigger', allow_none=True)
    obsNotes = List(Nested(ObsNoteSchema), attribute='obsNotes')
    exps = List(Nested(ExposureSchema), attribute='exps')
    telescopes = List(
        Nested(TelescopeBaseLimitedSchema), attribute='telescopes')


class DataProcessingJobLimitedSchema(SkynetSchema):
    id = Integer(attribute='id')


class DataProcessingFileSchema(SkynetSchema):
    id = Integer(attribute='id')
    jobId = Integer(attribute='jobId')
    type = String(attribute='type')
    state = String(attribute='state')

    job = Nested(DataProcessingJobLimitedSchema)

    dirPath = String(attribute='dirPath', dump_only=True)
    baseFilename = String(attribute='baseFilename', dump_only=True)
    filename = String(attribute='filename', dump_only=True)
    filePath = String(attribute='filerPath', dump_only=True)


class DataProcessingJobSchema(DataProcessingJobLimitedSchema):
    userId = Integer(attribute='userId')
    state = String(attribute='state')
    type = String(attribute='type')
    timeInUtc = DateTime(attribute='timeInUtc')
    timeFinishedUtc = DateTime(attribute='timeFinishedUtc')
    errorMsg = String(attribute='errorMsg')

    dirPath = String(attribute='dirPath', dump_only=True)

    files = List(Nested(DataProcessingFileSchema))


class RadioCartographerJobSchema(DataProcessingJobSchema):
    radioObsId = Integer(attribute='radioObsId')
    isRaw = Boolean(attribute='isRaw', allow_none=True)
    isLarge = Boolean(attribute='isLarge', allow_none=True)
    bgScale = Float(attribute='bgScale', allow_none=True)
    timeShift = Integer(attribute='timeShift', allow_none=True)
    rfiScale = Float(attribute='rfiScale', allow_none=True)
    weightScale = Float(attribute='weightScale', allow_none=True)
    noiseLevelPrior = Integer(attribute='noiseLevelPrior', allow_none=True)
    channel = String(attribute='channel')
    calibration = String(attribute='calibration', allow_none=True)
    mappingCoordinates = String(attribute='mappingCoordinates', allow_none=True)
    trimSize = Float(attribute='trimSize', allow_none=True)
    photometry = Integer(attribute='photometry', allow_none=True)
    aperatureRadius = Float(attribute='aperatureRadius')
    annulusRadius = Float(attribute='annulusRadius')
    sourceFindingMethod = String(attribute='sourceFindingMethod')
    coordinates = String(attribute='coordinates')

    radioObs = Nested(RadioObsLimitedSchema)


class RadioObsSchema(RadioObsLimitedSchema):
    priority = Integer(attribute='priority')
    state = String(attribute='state')
    dmsState = String(attribute='dmsState')
    dmsTransferAttempts = Integer(attribute='dmsTransferAttempts')
    timeIn = DateTime(attribute='timeIn')
    nChannel = Integer(attribute='nChannel')
    centerFrequency = Float(attribute='centerFrequency')
    integrationTime = Float(attribute='integrationTime')
    recordType = String(attribute='recordType')
    # autoPeak = Boolean(attribute='autoPeak')
    # autoFocus = Boolean(attribute='autoFocus')
    coordType = String(attribute='coordType')
    centerV = FloatSex(attribute='centerV', allow_none=True)
    centerH = FloatSex(attribute='centerH', allow_none=True)
    centerObjectName = String(attribute='centerObjectName', allow_none=True)
    centerObjectType = String(attribute='centerObjectType', allow_none=True)
    startAfter = DateTime(attribute='startAfter', allow_none=True)
    endBefore = DateTime(attribute='endBefore', allow_none=True)
    minEl = FloatSex(attribute='minEl', allow_none=True)
    timeCompleted = DateTime(attribute='timeCompleted', allow_none=True)
    minSolarSep = FloatSex(attribute='minSolarSep')
    # calInterval = Float(attribute='calInterval')
    remoteFilePath = String(attribute='remoteFilePath', allow_none=True)
    filePathNew = String(attribute='filePathNew', allow_none=True)
    numFiles = Integer(attribute='numFiles', allow_none=True)
    secondaryFrequency = Float(attribute='secondaryFrequency', allow_none=True)

    userId = Integer(attribute='userId', allow_none=True)
    groupId = Integer(attribute='groupId', allow_none=True)
    collabId = Integer(attribute='collabId', allow_none=True)
    receiverId = Integer(attribute='receiverId')
    radioFilterId = Integer(attribute='radioFilterId')
    timeAccountId = Integer(attribute='timeAccountId', allow_none=True)
    accountType = String(attribute='accountType', allow_none=True)
    teleOwnerId = Integer(attribute='teleOwnerId', allow_none=True)
    isToo = Boolean(attribute='isToo')
    currentAntennaId = Integer(attribute='currentAntennaId', allow_none=True)
    creditsCharged = Integer(attribute='creditsCharged', allow_none=True)

    leftJobId = Integer(attribute='leftJobId', allow_none=True)
    rightJobId = Integer(attribute='rightJobId', allow_none=True)
    compositeJobId = Integer(attribute='compositeJobId', allow_none=True)

    antennas = List(Nested(TelescopeBaseLimitedSchema), attribute='antennas')
    currentAntenna = Nested(
        TelescopeBaseLimitedSchema, attribute='currentAntenna', allow_none=True)
    radioFilter = Nested(
        RadioFilterLimitedSchema, attribute='radioFilter', allow_none=True)
    receiver = Nested(
        ReceiverLimitedSchema, attribute='receiver', allow_none=True)
    user = Nested(UserLimitedSchema, attribute='user', allow_none=True)
    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    collab = Nested(
        CollaborationLimitedSchema, attribute='collab', allow_none=True)
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    leftJob = Nested(
        RadioCartographerJobSchema, attribute='leftJob', allow_none=True)
    rightJob = Nested(
        RadioCartographerJobSchema, attribute='rightJob', allow_none=True)
    compositeJob = Nested(
        RadioCartographerJobSchema, attribute='compositeJob', allow_none=True)
    radioCartographerJobs = List(
        Nested(RadioCartographerJobSchema), attribute='radioCartographerJobs',
        allow_none=True)


class RadioObsTrackSchema(RadioObsSchema):
    duration = Float(attribute='duration')
    fixedOffsetH = FloatSex(attribute='fixedOffsetH')
    fixedOffsetV = FloatSex(attribute='fixedOffsetV')
    endOffsetH = FloatSex(attribute='endOffsetH')
    endOffsetV = FloatSex(attribute='endOffsetV')
    # cosv = Integer(attribute='cosv')
    repeat = Integer(attribute='repeat', default=0)
    obsId = Integer(attribute='obsId')


class RadioObsMapSchema(RadioObsSchema):
    hLength = FloatSex(attribute='hLength')
    vLength = FloatSex(attribute='vLength')
    delta = FloatSex(attribute='delta')
    duration = Float(attribute='duration')
    # cosv = Integer(attribute='cosv', default=1)
    # unidirectional = Integer(attribute='unidirectional', default=0)
    # start = Integer(attribute='start', default=0)
    # stop = Integer(attribute='stop', allow_none=True)
    mapType = String(attribute='mapType')
    obsId = Integer(attribute='obsId')

    pointsPerSweep = Float(attribute='pointsPerSweep', dump_only=True)
    gapBetweenSweeps = Float(attribute='gapBetweenSweeps', dump_only=True)
    vStar = Float(attribute='vStar', dump_only=True)
    gapAlongSweeps = Float(attribute='gapAlongSweeps', dump_only=True)
    slewSpeed = Float(attribute='slewSpeed', dump_only=True)
    numSweeps = Float(attribute='numSweeps', dump_only=True)


class RadioObsDaisySchema(RadioObsSchema):
    radius = Float(attribute='radius')
    radialPeriod = Float(attribute='radialPeriod')
    radialPhase = FloatSex(attribute='radialPhase', default=0)
    rotationPhase = FloatSex(attribute='rotationPhase', default=0)
    duration = Float(attribute='duration')
    # cosv = Integer(attribute='cosv', default=1)
    # calcDt = Float(attribute='calcDt', default=0.5)
    obsId = Integer(attribute='obsId')

    numPetals = Integer(attribute='numPetals', dump_only=True)


class RadioObsOnOffSchema(RadioObsSchema):
    duration = Float(attribute='duration')
    referenceOffsetH = FloatSex(attribute='referenceOffsetH')
    referenceOffsetV = FloatSex(attribute='referenceOffsetV')
    # cosv = Integer(attribute='cosv')
    reverse = Integer(attribute='reverse')
    repeat = Integer(attribute='repeat', default=0)
    obsId = Integer(attribute='obsId')


class TeleStatus2Schema(TeleStatus2LimitedSchema):
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)
    obs = Nested(ObservationLimitedSchema, attribute='obs', allow_none=True)
    exp = Nested(ExposureLimitedSchema, attribute='exp', allow_none=True)
    radioObs = Nested(
        RadioObsLimitedSchema, attribute='radioObs', allow_none=True)
    isValid = Boolean(attribute='isValid')


class SiteSchema(SiteLimitedSchema):
    sitePermissionGranted = Boolean(attribute='sitePermissionGranted')
    displayName = String(attribute='displayName', allow_none=True)
    locationName = String(attribute='locationName', allow_none=True)
    elevation = Integer(attribute='elevation', allow_none=True)
    skyBrightness = Float(attribute='skyBrightness', allow_none=True)
    avgSeeing = Float(attribute='avgSeeing', allow_none=True)
    teleSiteName = String(attribute='teleSiteName', allow_none=True)
    webcamUrl = String(attribute='webcamUrl', allow_none=True)
    fractionTimeObservable = Float(attribute='fractionTimeObservable')
    thumbnail = String(attribute='thumbnail', allow_none=True)
    wundergroundStationId = String(
        attribute='wundergroundStationId', allow_none=True)
    clearSkyClockId = String(attribute='clearSkyClockId', allow_none=True)
    description = String(attribute='description')

    telescopes = List(
        Nested(TelescopeBaseLimitedSchema), attribute='telescopes')
    antennas = List(
        Nested(TelescopeBaseLimitedSchema), attribute='antennas')
    admins = List(Nested(UserLimitedSchema), attribute='admins')


class FilterTransmissionSchema(SkynetSchema):
    __id_field__ = 'filterId'
    __name_field__ = None

    filterId = Integer(attribute='filterId')
    nu = Float(attribute='nu')
    transmission = Float(attribute='transmission')


class StandardFilterSchema(FilterSchema):
    zeroPoint = Float(attribute='zeroPoint', allow_none=True)
    flatOrder = Integer(attribute='flatOrder', allow_none=True)
    htmlHexColor = String(attribute='htmlHexColor', allow_none=True)
    transmissionData = List(
        Nested(FilterTransmissionSchema), attribute='transmissionData')


class GenericFilterSchema(FilterSchema):
    filters = List(Nested(StandardFilterSchema), attribute='filters')


class PriorityAccessScheduleSchema(PriorityAccessScheduleLimitedSchema):
    startDate = DateTime(attribute='startDate')
    stopDate = DateTime(attribute='stopDate', allow_none=True)
    userId = Integer(attribute='userId', allow_none=True)
    groupId = Integer(attribute='groupId', allow_none=True)
    collabId = Integer(attribute='collabId', allow_none=True)
    obsId = Integer(attribute='obsId', allow_none=True)
    order = Integer(attribute='order')
    teleOwnerId = Integer(attribute='teleOwnerId')
    teleId = Integer(attribute='teleId', allow_none=True)
    isToo = Boolean(attribute='isToo')

    collab = Nested(
        CollaborationLimitedSchema, attribute='collab', allow_none=True)
    group = Nested(GroupLimitedSchema, attribute='group', allow_none=True)
    user = Nested(UserLimitedSchema, attribute='user', allow_none=True)
    observation = Nested(
        ObservationLimitedSchema, attribute='observation', allow_none=True)
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)

    isActive = Boolean(attribute='isActive')


class TeleOwnerGuestShareSchema(TeleOwnerGuestShareLimitedSchema):
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)

    shareUsed = Float(attribute='shareUsed', dump_only=True)
    priority = Float(attribute='priority', dump_only=True)


class TeleOwnerShareSchema(TeleOwnerShareLimitedSchema):
    teleOwner = Nested(
        TeleOwnerLimitedSchema, attribute='teleOwner', allow_none=True)
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)

    contesting = Integer(attribute='contesting', dump_only=True)
    priorityAccessSchedules = List(
        Nested(PriorityAccessScheduleLimitedSchema),
        attribute='priorityAccessSchedules', dump_only=True)
    tooAccessSchedules = List(
        Nested(PriorityAccessScheduleLimitedSchema),
        attribute='tooAccessSchedules', dump_only=True)
    shareUsed = Float(attribute='shareUsed', dump_only=True)
    priority = Float(attribute='priority', dump_only=True)
    ownerShareUsed = Float(attribute='ownerShareUsed', dump_only=True)
    guestShareUsed = Float(attribute='guestShareUsed', dump_only=True)
    publicShareUsed = Float(attribute='publicShareUsed', dump_only=True)
    # noinspection PyTypeChecker
    accountTypeOrder = List(
        String, attribute='accountTypeOrder', dump_only=True)


class DarkCalibrationConfigLimitedSchema(SkynetSchema):
    __id_field__ = 'darkConfigId'
    __name_field__ = None

    darkConfigId = Integer(attribute='darkConfigId')
    teleId = Integer(attribute='teleId')
    numExps = Integer(attribute='numExps')
    expLength = Float(attribute='expLength')


class DarkCalibrationConfigSchema(DarkCalibrationConfigLimitedSchema):
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)


class FlatCalibrationConfigLimitedSchema(SkynetSchema):
    __id_field__ = 'flatConfigId'
    __name_field__ = None

    flatConfigId = Integer(attribute='flatConfigId')
    teleId = Integer(attribute='teleId')
    numExps = Integer(attribute='numExps')
    filterId = Integer(attribute='filterId')


class FlatCalibrationConfigSchema(FlatCalibrationConfigLimitedSchema):
    filter = Nested(FilterLimitedSchema, attribute='filter', allow_none=True)
    telescope = Nested(
        TelescopeBaseLimitedSchema, attribute='telescope', allow_none=True)


class TelescopeBaseSchema(TelescopeBaseLimitedSchema):
    teleType = String(attribute='teleType')
    ipAddress = String(attribute='ipAddress')
    siteId = Integer(attribute='siteId')
    thumbnail = String(attribute='thumbnail', allow_none=True)
    webcamUrl = String(attribute='webcamUrl', allow_none=True)
    description = String(attribute='description', allow_none=True)
    primaryImage = String(attribute='primaryImage', allow_none=True)
    isAvailable = Boolean(attribute='isAvailable')
    lastTeleOwnerUsageReset = Date(
        attribute='lastTeleOwnerUsageReset', allow_none=True)
    minEl = Float(attribute='minEl')
    maxSlewRateAxis1 = Float(attribute='maxSlewRateAxis1')
    maxSlewRateAxis2 = Float(attribute='maxSlewRateAxis2')

    localHorizonData = List(
        Nested(LocalHorizonSchema), attribute='localHorizonData')
    ownerShares = List(
        Nested(TeleOwnerShareLimitedSchema), attribute='ownerShares')
    guestShares = List(
        Nested(TeleOwnerGuestShareLimitedSchema), attribute='guestShares')
    status = Nested(
        TeleStatus2LimitedSchema, attribute='status', allow_none=True)
    admins = List(Nested(UserLimitedSchema), attribute='admins', dump_only=True)


class TelescopeLimitedSchema(TelescopeBaseLimitedSchema):
    aperture = Float(attribute='aperture', allow_none=True)
    focalLength = Float(attribute='focalLength', allow_none=True)
    pixScale = Float(attribute='pixScale')
    pixelScale = Float(attribute='pixelScale', dump_only=True)
    fov = List(Float(), attribute='fov', dump_only=True)
    filters = List(Nested(FilterLimitedSchema), attribute='filters')
    grbFilters = List(Nested(FilterLimitedSchema), attribute='grbFilters')
    genericFilters = List(
        Nested(FilterLimitedSchema), attribute='genericFilters')
    grbFilterNames = String(attribute='grbFilterNames')
    efficiency = Float(attribute='efficiency')
    flatEfficiency = Float(attribute='flatEfficiency')
    maxTrackingDuration = Float(attribute='maxTrackingDuration')


class TelescopeSchema(TelescopeBaseSchema, TelescopeLimitedSchema):
    tcsPort = Integer(attribute='tcsPort')
    settleTime = Float(attribute='settleTime')
    mountType = String(attribute='mountType')
    meridianSec = Integer(attribute='meridianSec')
    owner = String(attribute='owner')
    flatTarget = Float(attribute='flatTarget')
    flatTolerance = Float(attribute='flatTolerance')
    fileBlockSize = Integer(attribute='fileBlockSize')
    maxMbs = Float(attribute='maxMbs')
    flatMaxTime = Float(attribute='flatMaxTime')
    flatMinTime = Float(attribute='flatMinTime')
    darkElev = Float(attribute='darkElev')
    flatElev = Float(attribute='flatElev')
    lightElev = Float(attribute='lightElev')
    dmsPort = Integer(attribute='dmsPort')
    extraText = String(attribute='extraText')
    dmsIP = String(attribute='dmsIP', allow_none=True)
    avgWaitTime = Float(attribute='avgWaitTime', allow_none=True)
    teleGroupId = Integer(attribute='teleGroupId', allow_none=True)
    rbiTimeConstant = Float(attribute='rbiTimeConstant', allow_none=True)
    rbiAvgBkg = Integer(attribute='rbiAvgBkg', allow_none=True)
    allowBrightTargets = Boolean(attribute='allowBrightTargets')
    scheduleNextFlatObsOn = DateTime(
        attribute='scheduleNextFlatObsOn', allow_none=True)
    scheduleNextDarkObsOn = DateTime(
        attribute='scheduleNextDarkObsOn', allow_none=True)
    allowObsLinking = Boolean(attribute='allowObsLinking')
    initiateDmsConnection = Boolean(attribute='initiateDmsConnection')
    enforceMaxExpLength = Boolean(attribute='enforceMaxExpLength')
    allowCustomBinning = Boolean(attribute='allowCustomBinning')
    defaultBinning = Integer(attribute='defaultBinning')

    dmsStatus = Nested(DmsStatusSchema, attribute='dmsStatus', allow_none=True)
    darkCalibrationConfigs = List(
        Nested(DarkCalibrationConfigLimitedSchema),
        attribute='darkCalibrationConfigs')
    flatCalibrationConfigs = List(
        Nested(FlatCalibrationConfigLimitedSchema),
        attribute='flatCalibrationConfigs')
    # TODO: cameraHistory = relationship(u'CameraHistory', backref='telescope')

    # TODO: camera = Nested(CameraLimitedSchema,
    #     attribute='camera', allow_none=True)

    uniqueId = String(attribute='uniqueId', dump_only=True)


class ReceiverSchema(ReceiverLimitedSchema):
    bandStartFreq = Float(attribute='bandStartFreq')
    bandStopFreq = Float(attribute='bandStopFreq')
    hasTunableLo = Boolean(attribute='hasTunableLo')
    discreteFrequencyOffset = Float(attribute='discreteFrequencyOffset')
    lowResBandwidth = Float(attribute='lowResBandwidth')
    lowResUsesFilters = Boolean(attribute='lowResUsesFilters')
    lowResDualBand = Boolean(attribute='lowResDualBand')
    lowResCenterFreqDefault = Float(
        attribute='lowResCenterFreqDefault', allow_none=True)
    lowResSecondaryFreqDefault = Float(
        attribute='lowResSecondaryFreqDefault', allow_none=True)
    highResBandwidth = Float(attribute='highResBandwidth')
    highResUsesFilters = Boolean(attribute='highResUsesFilters')
    highResDualBand = Boolean(attribute='highResDualBand')
    highResCenterFreqDefault = Float(
        attribute='highResCenterFreqDefault', allow_none=True)
    highResSecondaryFreqDefault = Float(
        attribute='highResSecondaryFreqDefault', allow_none=True)
    filters = List(Nested(RadioFilterSchema), attribute='filters')


class ReceiverHistoryLimitedSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    teleId = Integer(attribute='teleId')
    receiverId = Integer(attribute='receiverId')
    startDate = DateTime(attribute='startDate')


class ReceiverHistorySchema(ReceiverHistoryLimitedSchema):
    receiver = Nested(
        ReceiverLimitedSchema, attribute='receiver', allow_none=True)
    antenna = Nested(
        TelescopeBaseLimitedSchema, attribute='antenna', allow_none=True)


class AntennaLimitedSchema(TelescopeBaseLimitedSchema):
    diameter = Float(attribute='diameter')


class AntennaSchema(TelescopeBaseSchema, AntennaLimitedSchema):
    slewOverhead = Float(attribute='slewOverhead')
    dmsIp = String(attribute='dmsIp', allow_none=True)

    receiverHistory = List(
        Nested(ReceiverHistoryLimitedSchema), attribute='receiverHistory')

    receiver = Nested(ReceiverSchema, attribute='receiver', allow_none=True)


class TeleOwnerSchema(TeleOwnerLimitedSchema):
    priorityAccessSchedules = List(
        Nested(PriorityAccessScheduleSchema),
        attribute='priorityAccessSchedules')
    shares = List(Nested(TeleOwnerShareLimitedSchema), attribute='shares')
    admins = List(
        Nested(UserLimitedSchema), attribute='admins', dump_only=True)


class TimeAccountTxnSchema(SkynetSchema):
    __name_field__ = None

    id = Integer(attribute='id')
    timeAccountId = Integer(attribute='timeAccountId')
    amount = Integer(attribute='amount')
    comment = String(attribute='comment', allow_none=True)
    time = DateTime(attribute='time')
    obsId = Integer(attribute='obsId', allow_none=True)
    expId = Integer(attribute='expId', allow_none=True)
    partnerTxnId = Integer(attribute='partnerTxnId', allow_none=True)
    partnerTimeAccountId = Integer(attribute='partnerTimeAccountId')


class GroupSchema(GroupLimitedSchema):
    priority = Integer(attribute='priority')
    owningCollabId = Integer(attribute='owningCollabId', allow_none=True)
    owningCollab = Nested(
        CollaborationLimitedSchema, attribute='owningCollab', allow_none=True)
    users = List(Nested(UserLimitedSchema), attribute='users')
    ownedUsers = List(Nested(UserLimitedSchema), attribute='ownedUsers')
    timeAccounts = List(
        Nested(TimeAccountLimitedSchema), attribute='timeAccounts')

    collabs = List(Nested(CollaborationLimitedSchema), attribute='collabs')
    pendingGroupMemberships = List(
        Nested(PendingGroupMembershipLimitedSchema),
        attribute='pendingGroupMemberships')
    pendingCollabMemberships = List(
        Nested(PendingCollabMembershipLimitedSchema),
        attribute='pendingCollabMemberships')

    admins = List(Nested(UserLimitedSchema), attribute='admins', dump_only=True)


class CollaborationSchema(CollaborationLimitedSchema):
    groups = List(Nested(GroupLimitedSchema), attribute='groups')
    timeAccounts = List(
        Nested(TimeAccountLimitedSchema), attribute='timeAccounts')

    ownedGroups = List(Nested(GroupLimitedSchema), attribute='ownedGroups')
    pendingCollabMemberships = List(
        Nested(PendingCollabMembershipLimitedSchema),
        attribute='pendingCollabMemberships')

    admins = List(Nested(UserLimitedSchema), attribute='admins', dump_only=True)


class JobSchema(SkynetSchema):
    id = Integer(attribute='id')
    userId = Integer(attribute='userId')
    apiVersion = String(attribute='apiVersion')
    status = String(attribute='status')
    timeCreated = DateTime(attribute='timeCreated')
    timeCompleted = DateTime(attribute='timeCompleted', allow_none=True)
    type = String(attribute='type')
    args = String(attribute='args')
    result = String(attribute='result', allow_none=True)

    user = Nested(UserLimitedSchema, attribute='user')


__all__ = ['is_scalar_field'] + [
    _name for _name, _val in globals().items() if isinstance(_val, type) and (
        issubclass(_val, fields.Field) or
        issubclass(_val, Schema) and _val is not Schema)]
