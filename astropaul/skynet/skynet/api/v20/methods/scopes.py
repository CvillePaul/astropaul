"""
SkyNet API procedural Python interface: telescope methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get']


def query(site=None, sort_by=None, offset=None, limit=None, include=None,
          exclude=None, api_key=None, **kwargs):
    """
    Return a list of telescopes matching certain criteria; the user must be
    a Skynet or telescope admin or have a time account that provides access to
    some telescopes

    :param int | str | list site: return only telescopes located at the given
        site(s); must be a single integer site ID, a string containing a
        comma-separated list of site IDs or names, or a list of IDs or names
    :param str sort_by: sort telescopes by the given column(s); should be
        a comma-separated list of TelescopeBase, Telescope, or Antenna column
        names, with "-" before the column name indicating reversed order
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
    :param kwargs: any other TelescopeBase field=value pairs to filter by

    :return: list of serialized :class:`skynet.api.serializers.TelescopeSchema`
        or :class:`skynet.api.serializers.AntennaSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        site=site, sort_by=sort_by, offset=offset, limit=limit, include=include,
        exclude=exclude))
    return api_call('scopes', 'get', kwargs, api_key=api_key)


def get(scope, api_key=None):
    """
    Return the given telescope profile; only telescope admins can access the
    full profiles; users with time accounts that provide access to the given
    telescope can view limited profiles

    :param int | str scope: telescope ID or name
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.orm.Telescope` or
        :class:`skynet.api.orm.Antenna` instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('scopes/{}'.format(scope), api_key=api_key)
