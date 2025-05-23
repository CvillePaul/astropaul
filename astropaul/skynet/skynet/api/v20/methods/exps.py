"""
SkyNet API procedural Python interface: exposure methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'add', 'update']


def query(obs=None, user=None, group=None, scope=None, wcs_present=None,
          after=None, before=None, sort_by=None, offset=None, limit=None,
          include=None, exclude=None, api_key=None, **kwargs):
    """
    Return a list of exposures matching certain criteria

    :param int | str | list obs: only return exposures for specific observation
        ID(s) or name(s); must be a single integer observation ID, a string
        containing a comma-separated list of observation IDs or names, or a list
        of IDs or names
    :param int | str | list user: for SkyNet admins and group admins only:
        request exposures submitted by specific user(s); must be a single
        integer user ID, a comma-separated list of user IDs or usernames, or a
        list of user IDs or usernames; default: return user's own exposures
    :param int | str | list group: for SkyNet admins and group admins only:
        request exposures submitted by all users of specific group(s); must be a
        single integer group ID, a string containing a comma-separated list of
        group IDs or names, or a list of group IDs or names
    :param int | str | list scope: only return exposures taken with specific
        telescope(s); must be a single integer telescope ID, a string containing
        a comma-separated list of telescope IDs/names, or a list of IDs or names
    :param str | int wcs_present: limit to exposures with (1) or without (0) a
        WCS solution; implies wcs_state="completed"
    :param DateTime | str after: request exposures taken after the given
        date/time
    :param DateTime | str before: request exposures taken before the given
        date/time
    :param str sort_by: sort exposures by the given column(s); should be
        a comma-separated list of Exposure column names, with "-" before
        the column name indicating reversed order
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
    :param kwargs: any other Exposure field=value pairs to filter by

    Note. Setting either `obs`, `user`, or `group` allows the user to get
    observations submitted by other users, provided the user has privileges to
    view others' observations; these include the collaboration, group, and
    telescope admins of the observation's time account and SkyNet admins. If
    neither of those parameters are set, the user gets only their own
    exposures, possibly restricted by other parameters. To get _all_ exposures
    for the given telescope, admin needs to set one of those parameters to match
    all exposures, e.g. scope='Morehead', user=''.

    :return: list of serialized :class:`skynet.api.serializers.ExposureSchema`
        instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        obs=obs, user=user, group=group, scope=scope, wcs_present=wcs_present,
        after=after, before=before, sort_by=sort_by, offset=offset,
        include=include, exclude=exclude, limit=limit))
    return api_call('exps', 'get', kwargs, api_key=api_key)


def get(exp, api_key=None):
    """
    Return exposure with the given ID; the calling user must be either the
    observation submitter, the observation's time account collaboration, group,
    or telescope admin, or a SkyNet admin

    :param int | str exp: exposure ID
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.ExposureSchema` instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('exps/{}'.format(exp), api_key=api_key)


def add(obs, api_key=None, **kwargs):
    """
    Add an exposure to an existing observation

    :param int | str obs: ID of the observation to update
    :param str api_key: optional API access key
    :param kwargs: any skynet.api.serializers.ExposureSchema fields

    :return: serialized skynet.api.serializers.ExposureSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    kwargs = dict(kwargs)
    kwargs['obs'] = obs
    return api_call('exps', 'post', kwargs, api_key=api_key)


def update(exp, api_key=None, **kwargs):
    """
    Modify the existing exposure

    :param int | str exp: ID of the exposure to update
    :param str api_key: optional API access key
    :param kwargs: any skynet.api.serializers.ExposureSchema fields to set

    :return: serialized skynet.api.serializers.ExposureSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('exps/{}'.format(exp), 'put', kwargs, api_key=api_key)
