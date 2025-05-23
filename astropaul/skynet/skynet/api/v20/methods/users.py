"""
SkyNet API procedural Python interface: user methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'create', 'update', 'delete']


def query(group=None, sort_by=None, offset=None, limit=None, include=None,
          exclude=None, api_key=None, **kwargs):
    """
    Return a list of users matching certain criteria; except for SkyNet admins,
    only users owned by or belonging to the groups administered by the calling
    user and the calling user's ID are returned

    :param int | str | list group: request users belonging to the given
        group(s); must be a single integer group ID, a string containing a
        comma-separated list of group IDs or names, or a list of IDs or names;
        the user must be a group or an owning collab admin; SkyNet admins must
        set `group`='' to get all SkyNet users
    :param str sort_by: sort users by the given column(s); should be
        a comma-separated list of User column names, with "-" before the column
        name indicating reversed order
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
    :param kwargs: any other User field=value pairs to filter by

    :return: list of serialized
        :class:`skynet.api.serializers.UserSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        group=group, sort_by=sort_by, offset=offset, limit=limit,
        include=include, exclude=exclude))
    return api_call('users', 'get', kwargs, api_key=api_key)


def get(user=None, api_key=None):
    """
    Return user with the given ID or name; only SkyNet admins and users
    themselves can get the full user profiles; owning group admins can get
    partial user profiles with sensitive information hidden; group admins can
    get IDs and usernames of group members

    :param str | int user: user ID or name; defaults to the calling user
    :param str api_key: optional API access key

    :return: serialized skynet.api.serializers.UserSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('users/{}'.format(user) if user is not None else 'users',
                    api_key=api_key)


def create(username, email, api_key=None, **kwargs):
    """
    Create a new user

    :param str username: SkyNet username
    :param str email: user's email address
    :param str api_key: optional API access key
    :param kwargs:
        owningGroup: ID, name, or API ORM object of the owning group; if given
            and non-empty, the calling user must be that group's admin with a
            CREATE_MEMBERS role or a SkyNet admin; users with no owning group
            can be created by SkyNet admins only

    :return: serialized skynet.api.serializers.UserSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    kwargs = dict(kwargs)
    kwargs['username'] = username
    kwargs['email'] = email
    return api_call('users', 'post', kwargs, api_key=api_key)


def update(user=None, api_key=None, **kwargs):
    """
    Modify the existing user data

    :param int | str user: ID or username of the user to update; group admins
        with a MANAGE_MEMBERS role for the user's owning group can modify any
        non-personal scalar fields; when updating the calling user's own
        profile, non-admin user can modify only personal data and preferences;
        if unspecified, update the calling user's profile
    :param str api_key: optional API access key
    :param kwargs: any skynet.api.serializers.UserSchema fields to set

    :return: serialized skynet.api.serializers.UserSchema instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('users/{}'.format(user) if user is not None else 'users',
                    'put', kwargs, api_key=api_key)


def delete(user, api_key=None):
    """
    Permanently delete a user

    Currently only a user can be deleted who has no associated time accounts
    and observations. Only SkyNet admins and admins of all user's groups with
    MANAGE_MEMBERS role can do that.

    :param int | str user: ID or username
    :param str api_key: optional API access key

    :return: None
    """
    return api_call('users/{}'.format(user), 'delete', api_key=api_key)
