"""
SkyNet API procedural Python interface: filter methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get']


def query(api_key=None, **kwargs):
    """
    Return a list of filters matching certain criteria

    :param str api_key: optional API access key
    :param kwargs: any Filter field=value pairs to filter by, plus
        sort_by: sort filters by the given column(s); should be
            a comma-separated list of Filter column names, with "-" before
            the column name indicating reversed order
        offset: only return items starting from the given index;
            default: start from 0
        limit: only return the given number of items at maximum;
            default: return no more than 1000 items
        include: list of fields to return; default: return only "id" if
            `exclude` is not set, otherwise return all fields except those
            listed in `exclude`; the sole "include=*" returns all fields
        exclude: list of fields to exclude from serialization; default: return
            all fields listed in `include`; `exclude` takes precedence if
            a field is listed both in `include` and `exclude`

    :return: list of serialized
        :class:`skynet.api.serializers.GenericFilterSchema` or
        :class:`skynet.api.serializers.StandardFilterSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    return api_call('filters', 'get', kwargs, api_key=api_key)


def get(flt, api_key=None):
    """
    Return filter with the given ID or name

    :param int | str flt: filter ID or name
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.GenericFilterSchema` or
        `StandardFilterSchema` instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('filters/{}'.format(flt), api_key=api_key)
