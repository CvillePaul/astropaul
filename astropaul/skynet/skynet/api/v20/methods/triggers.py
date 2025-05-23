"""
SkyNet API procedural Python interface: trigger methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'add', 'update']


def query(user=None, sort_by=None, offset=None, limit=None, include=None,
          exclude=None, api_key=None, **kwargs):
    """
    Return a list of campaign triggers matching certain criteria

    :param int | str | list user: for SkyNet admins and group admins only:
        request campaigns submitted by specific user(s); must be a single
        integer user ID, a string containing a comma-separated list of user IDs
        or usernames, or a list of IDs or usernames; default: return user's own
        campaigns
    :param str sort_by: sort triggers by the given column(s); should be
        a comma-separated list of CampaignTrigger column names, with "-"
        before the column name indicating reversed order
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
    :param kwargs: any other Campaign field=value pairs to filter by

    Note. Setting `user` allows the user to get campaigns
    submitted by other users, provided the user has privileges to view others'
    campaigns; these include SkyNet admins. If neither of those parameters
    are set, the user gets only their own campaigns, possibly restricted
    by other parameters.

    :return: list of serialized
        :class:`skynet.api.serializers.CampiagnSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        user=user, sort_by=sort_by, offset=offset, limit=limit, include=include,
        exclude=exclude))
    return api_call('triggers', 'get', kwargs, api_key=api_key)


def get(trigger, api_key=None):
    """
    Return trigger with the given ID or name; the calling user must
    be either the trigger submitter or a SkyNet admin

    :param int | str trigger: trigger ID or name
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.CampaignTriggerSchema`
        instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('triggers/{}'.format(trigger), api_key=api_key)


def add(api_key=None, **kwargs):
    """
    Submit a trigger

    :param str api_key: optional API access key
    :param kwargs: any :class:`skynet.api.serializers.CampiagnTriggerSchema`
        fields

    :return: serialized :class:`skynet.api.serializers.CampaignTriggerSchema`
        instance for the new trigger
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('triggers', 'post', kwargs, api_key=api_key)


def update(trigger, api_key=None, **kwargs):
    """
    Modify the existing trigger

    :param int | str trigger: ID of the campaign to update
    :param str api_key: optional API access key
    :param kwargs: any :class:`skynet.api.serializers.CampaignTriggerSchema`
        fields to set

    :return: serialized :class:`skynet.api.serializers.CampaignTriggerSchema`
        instance
    :rtype: :class:`skynet.api.client.AttrDict`
    """
    return api_call(
        'triggers/{}'.format(trigger), 'put', kwargs, api_key=api_key)
