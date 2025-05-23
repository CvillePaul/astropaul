"""
SkyNet API procedural Python interface: collab methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'create', 'update', 'addgroup', 'delgroup', 'delete']


def query(group=None, sort_by=None, offset=None, limit=None, include=None,
          exclude=None, api_key=None, **kwargs):
    """
    Return a list of collaborations matching certain criteria; only collabs that
    the user has an admin role of and the user's own collabs are returned

    :param int | str | list group: request collabs of the given group(s); must
        be a single integer group ID, a string containing a comma-separated list
        of group IDs or names, or a list of IDs or names; default: return user's
        own collabs; SkyNet admins must set `group`='' to get all collabs
    :param str sort_by: sort collaborations by the given column(s); should be
        a comma-separated list of Collaboration column names, with "-" before
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
    :param kwargs: any other Collaboration field=value pairs to filter by

    :return: list of serialized
        :class:`skynet.api.serializers.CollaborationSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        group=group, sort_by=sort_by, offset=offset, limit=limit,
        include=include, exclude=exclude))
    return api_call('collabs', 'get', kwargs, api_key=api_key)


def get(collab, api_key=None):
    """
    Return collaboration with the given ID or name; only SkyNet and collab
    admins can get the full collab profiles; non-admin users can get a limited
    profile of a collab they belong to via any of their groups

    :param str | int collab: collab ID or name
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.CollaborationSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('collabs/{}'.format(collab), api_key=api_key)


def create(name, api_key=None):
    """
    Create a new empty collaboration

    :param str name: collaboration name
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.CollaborationSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('collabs', 'post', dict(name=name), api_key=api_key)


def update(collab, api_key=None, **kwargs):
    """
    Modify the existing collaboration data

    :param int | str collab: ID or name of the collaboration to update; the
        calling user must have an ADMIN role for the collab
    :param str api_key: optional API access key
    :param kwargs: any skynet.api.serializers.CollaborationSchema fields to set

    :return: serialized skynet.api.serializers.CollaborationSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('collabs/{}'.format(collab), 'put', kwargs, api_key=api_key)


def addgroup(collab, group, api_key=None):
    """
    Add group to collaboration

    :param int | str collab: ID or name of the collaboration; the calling user
        must have a MANAGE_MEMBERS role for the collab
    :param int | str group: ID or name of the group to add
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.CollaborationSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call(
        'collabs/{}/groups'.format(collab), 'post', dict(group=group),
        api_key=api_key)


def delgroup(collab, group, api_key=None):
    """
    Remove group from collaboration

    :param int | str collab: ID or name of the collaboration; the calling user
        must have a MANAGE_MEMBERS role for the collab
    :param int | str group: ID or name of the group to remove
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.CollaborationSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call(
        'collabs/{}/groups/{}'.format(collab, group), 'delete', api_key=api_key)


def delete(collab, api_key=None):
    """
    Permanently delete collaboration

    :param int | str collab: ID or name of the collaboration; the calling user
        must be a SkyNet admin
    :param str api_key: optional API access key

    :rtype: None
    """
    return api_call('collabs/{}'.format(collab), 'delete', api_key=api_key)
