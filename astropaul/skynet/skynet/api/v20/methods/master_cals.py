"""
SkyNet API procedural Python interface: master calibration image methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get']


def query(scope=None, obs=None, exp=None, delta=None, after=None, before=None,
          sort_by=None, offset=None, limit=None, include=None, exclude=None,
          api_key=None, **kwargs):
    """
    Return a list of master cals matching certain criteria

    :param str | int | list scope: only return masters for specific
        telescope(s); must be a single integer telescope ID, a string containing
        a comma-separated list of telescope IDs or names, or a list of IDs or
        names; users other than SkyNet admins may only get masters for
        telescopes they are admins of or those they have access to via any of
        their active time accounts; by default, SkyNet admins get masters for
        the telescopes they have a time account for unless they explicitly set
        `scope`
    :param str | int | list obs: only return masters for specific
        observation(s); must be a single integer obs ID, a string containing
        a comma-separated list of observation IDs or names, or a list of IDs or
        names
    :param str | int | list exp: only return masters for specific exposure(s);
        must be a single integer obs ID, a string containing a comma-separated
        list of exposure IDs, or a list of IDs or names
    :param str | float delta: maximum separation in days between master cal and
        exposure epochs for `obs` and `exp`; default: 10
    :param str | datetime.datetime after: only return masters with startDate
        after the given date/time
    :param str | datetime.datetime before: only return masters with stopDate
        before the given date/time
    :param str sort_by: sort master cals by the given column(s); should be
        a comma-separated list of MasterCalibration column names, with "-"
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
    :param kwargs: any other MasterCalibration field=value pairs to filter by

    :return: list of serialized
        :class:`skynet.api.serializers.MasterCalibrationSchema` instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    kwargs = dict(kwargs)
    kwargs.update(dict(
        scope=scope, obs=obs, exp=exp, delta=delta, after=after, before=before,
        sort_by=sort_by, offset=offset, limit=limit, include=include,
        exclude=exclude))
    return api_call('master-cals', 'get', kwargs, api_key=api_key)


def get(id, api_key=None):
    """
    Return master cal with the given ID; the calling user must have access to
    the telescope by means of any of their time accounts or be a telescope admin

    :param int | str id: master cal ID
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.MasterCalibrationSchema`
        instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('master-cals/{}'.format(id), api_key=api_key)
