"""
SkyNet API procedural Python interface: group methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'create', 'update', 'adduser', 'deluser', 'delete']


def query(collab=None, user=None, sort_by=None, offset=None, limit=None,
          include=None, exclude=None, api_key=None, **kwargs):
    """
    Return a list of groups matching certain criteria; only groups that the user
    has an admin role for, groups belonging to collab administered by the user,
    and the user's own groups are returned

    :param int | str | list collab: request groups belonging to the given
        collab(s); must be a single integer collab ID, a string containing a
        comma-separated list of collab IDs or names, or a list of IDs or names;
        the user must be a SkyNet or collab admin
    :param int | str | list user: request groups of the given user(s); must be a
        single integer user ID, a string containing a comma-separated list of
        user IDs or usernames, or a list of IDs or usernames; each user must
        belong to a group administered by the calling user; SkyNet admins must
        set e.g. `user`='' to get all groups
    :param str sort_by: sort groups by the given column(s); should be
        a comma-separated list of Group column names, with "-" before
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
    :param kwargs: any other Group field=value pairs to filter by

    :return: list of serialized :class:`skynet.api.serializers.GroupSchema`
        instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        collab=collab, user=user, sort_by=sort_by, offset=offset, limit=limit,
        include=include, exclude=exclude))
    return api_call('groups', 'get', kwargs, api_key=api_key)


def get(group, api_key=None):
    """
    Return group with the given ID or name; only SkyNet and group admins can get
    the full collab profiles; non-admin group members and group's collab admins
    can get a limited group profile

    :param str | int group: group ID or name
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.GroupSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('groups/{}'.format(group), api_key=api_key)


def create(name, api_key=None, **kwargs):
    """
    Create a new group

    The user who created the group automatically joins it with an ADMIN role.

    :param str name: group name
    :param str api_key: optional API access key
    :param kwargs:
        owningCollab: ID or name of the owning collaboration; if given and
            non-empty, the calling user must be that collab's admin with a
            CREATE_MEMBERS role or a SkyNet admin; only SkyNet admins can create
            groups with no owning collaboration

    :return: serialized skynet.api.serializers.GroupSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    kwargs = dict(kwargs)
    kwargs['name'] = name
    return api_call('groups', 'post', kwargs, api_key=api_key)


def update(group, api_key=None, **kwargs):
    """
    Modify the existing group data

    :param int | str group: ID or name of the group to update; the calling user
        should have an ADMIN role for the group or for the group's owning collab
    :param str api_key: optional API access key
    :param kwargs: any skynet.api.serializers.GroupSchema fields to set

    :return: serialized skynet.api.serializers.GroupSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('groups/{}'.format(group), 'put', kwargs, api_key=api_key)


def adduser(group, user, api_key=None):
    """
    Add user to group

    :param int | str group: ID or name of the group; the calling user should
        have a MANAGE_MEMBERS role for the group
    :param int | str user: ID or name of the user to add
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.GroupSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call(
        'groups/{}/users'.format(group), 'post', dict(user=user),
        api_key=api_key)


def deluser(group, user, api_key=None):
    """
    Remove user from group

    :param int | str group: ID or name of the group; the calling user should
        have a MANAGE_MEMBERS role for the group
    :param int | str user: ID or name of the user to remove
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.GroupSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call(
        'groups/{}/users/{}'.format(group, user), 'delete', api_key=api_key)


def delete(group, api_key=None):
    """
    Permanently delete group

    :param int | str group: ID or name of the group; the calling user should
        have an ADMIN role for the group and a MANAGE_MEMBERS role for all
        collabs the group belongs to
    :param str api_key: optional API access key

    :rtype: None
    """
    return api_call('groups/{}'.format(group), 'delete', api_key=api_key)
