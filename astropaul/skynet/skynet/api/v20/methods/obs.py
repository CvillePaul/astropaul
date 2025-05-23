"""
SkyNet API procedural Python interface: observation methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'add', 'update']


def query(user=None, group=None, scope=None, after=None, before=None,
          sort_by=None, offset=None, limit=None, include=None, exclude=None,
          api_key=None, **kwargs):
    """
    Return a list of optical observations matching certain criteria

    :param int | str | list user: for SkyNet admins and group admins only:
        request observations submitted by specific user(s); must be a single
        integer user ID, a string containing a comma-separated list of user IDs
        or usernames, or a list of IDs or usernames; default: return user's own
        observations
    :param int | str | list group: for SkyNet admins and group admins only:
        request observations submitted by all users of specific group(s); must
        be a single integer group ID, a string containing a comma-separated list
        of group IDs or names, or a list of IDs or names; default: return user's
        own observations
    :param int | str | list scope: only return observations submitted to
        specific telescope(s); must be a single integer telescope ID, a string
        containing a comma-separated list of telescope IDs or names, or a list
        of IDs or names
    :param DateTime | str after: request observations submitted after the given
        date/time
    :param DateTime | str before: request observations submitted before the
        given date/time
    :param str sort_by: sort observations by the given column(s); should be
        a comma-separated list of Observation or Exposure column names
        (prepend "Exposure." if ambiguous), with "-" before the column name
        indicating reversed order
    :param int | str offset: only return items starting from the given index;
        default: start from 0
    :param str | int limit: only return the given number of items at maximum;
        default: return no more than 1000 items
    :param include: list of fields to return; default: return only "id" if
        `exclude` is not set, otherwise return all fields except those listed in
        `exclude`; the sole "include=*" returns all fields
    :param exclude: list of fields to exclude from serialization; default:
        return all fields listed in `include`; `exclude` takes precedence if
        a field is listed both in `include` and `exclude`
    :param str api_key: optional API access key
    :param kwargs: any other Observation field=value pairs to filter by

    Note. Setting either `user` or `group` allows the user to get observations
    submitted by other users, provided the user has privileges to view others'
    observations; these include the collaboration, group, and telescope admins
    of the observation's time account and SkyNet admins. If neither of those
    parameters are set, the user gets only their own observations, possibly
    restricted by other parameters. To get _all_ observations for the given
    telescope, admin needs to set one of those parameters to match all
    observations, e.g. scope='Morehead', user=''.

    :return: list of serialized
        :class:`skynet.api.serializers.ObservationSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        user=user, group=group, scope=scope, after=after, before=before,
        sort_by=sort_by, offset=offset, limit=limit, include=include,
        exclude=exclude))
    return api_call('obs', 'get', kwargs, api_key=api_key)


def get(obs, api_key=None):
    """
    Return optical observation with the given ID or name; the calling user must
    be either the observation submitter,the observation's time account
    collaboration, group, or telescope admin, or a SkyNet admin

    :param int | str obs: observation ID or name
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.ObservationSchema`
        instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('obs/{}'.format(obs), api_key=api_key)


def add(api_key=None, **kwargs):
    """
    Submit an observation

    :param str api_key: optional API access key
    :param kwargs: any :class:`skynet.api.serializers.ObservationSchema` fields,
        including scalar fields and such non-scalar fields as `telescopes`.
        Non-scalar fields must be given as dictionaries of the corresponding
        :mod:`skynet.api.serializers` schema fields.

    :return: serialized :class:`skynet.api.serializers.ObservationSchema`
        instance for the new observation
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('obs', 'post', kwargs, api_key=api_key)


def update(obs, api_key=None, **kwargs):
    """
    Modify the existing observation

    :param int | str obs: ID of the observation to update
    :param str api_key: optional API access key
    :param kwargs: any :class:`skynet.api.serializers.ObservationSchema` fields
        to set

    :return: serialized :class:`skynet.api.serializers.ObservationSchema`
        instance
    :rtype: :class:`skynet.api.client.AttrDict`
    """
    return api_call('obs/{}'.format(obs), 'put', kwargs, api_key=api_key)
