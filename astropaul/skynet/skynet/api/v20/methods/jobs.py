"""
SkyNet API procedural Python interface: job methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['query', 'get', 'get_status', 'get_result', 'add']


def query(api_key=None, **kwargs):
    """
    Return a list of jobs matching certain criteria; for anyone except Skynet
    admins, only the user's own jobs are returned

    :param str api_key: optional API access key
    :param kwargs: any Job field=value pairs to filter by, plus
        user: request the given user's jobs; must be a single integer user ID,
            a string containing a comma-separated list of user IDs or usernames,
            or a list of IDs or usernames; Skynet admins must set e.g. `user`=''
            to get jobs submitted by all users
        sort_by: sort jobs by the given column(s); should be
            a comma-separated list of Job column names, with "-" before
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

    :return: list of serialized :class:`skynet.api.serializers.JobSchema`
        instances
    :rtype: list[skynet.api.client.AttrDict]
    """
    return api_call('jobs', 'get', kwargs, api_key=api_key)


def get(job, api_key=None):
    """
    Return job with the given ID; the calling user must be a Skynet admin or
    the job's owner

    :param int job: job ID
    :param str api_key: optional API access key

    :return: serialized :class:`skynet.api.serializers.JobSchema` instance
    :rtype: skynet.api.client.AttrDict
    """
    return api_call('jobs/{}'.format(job), api_key=api_key)


def get_status(job, api_key=None):
    """
    Return the given job's status; the calling user must be a Skynet admin or
    the job's owner

    :param int job: job ID
    :param str api_key: optional API access key

    :return: job status ("pending", "running", "success", or "error")
    :rtype: str
    """
    return api_call('jobs/{}/status'.format(job), api_key=api_key)


def get_result(job, api_key=None):
    """
    Return the given job's result; the calling user must be a Skynet admin or
    the job's owner

    :param int job: job ID
    :param str api_key: optional API access key

    :return: job result if success, error message otherwise
    """
    return api_call('jobs/{}/result'.format(job), api_key=api_key)


def add(type, args, api_key=None):
    """
    Submit a job

    :param str type: job type
    :param str args: JSON-serialized job arguments
    :param str api_key: optional API access key

    :return: job ID
    :rtype: int
    """
    return api_call('jobs', 'post', dict(type=type, args=args), api_key=api_key)
