"""
SkyNet API: procedural Python interface core
"""

from __future__ import absolute_import, division, print_function

import sys
import os.path
import datetime
import warnings

import requests

from . import errors, __api_version__

if sys.version_info.major >= 3:
    from pickle import loads
    # noinspection PyCompatibility
    import configparser
else:
    # noinspection PyCompatibility,PyUnresolvedReferences
    from cPickle import loads
    # noinspection PyCompatibility,PyUnresolvedReferences
    import ConfigParser as configparser


__all__ = ['api_call', 'config', 'convert_json', 'AttrDict']


config = configparser.ConfigParser()
config.read(os.path.normpath(os.path.expanduser('~/.skynet/api.cfg')))
try:
    DEFAULT_SERVER = config.get('API', 'server')
except (configparser.NoSectionError, configparser.NoOptionError):
    DEFAULT_SERVER = 'production'

try:
    with open(os.path.normpath(os.path.expanduser('~/.skynet/API_KEY')),
              'r') as f:
        API_KEY = f.read().strip()
except (IOError, OSError):
    API_KEY = None

BASE_URL = {
    'local': 'http://127.0.0.1:8000/api{}'.format(__api_version__),
    'local_https': 'https://127.0.0.1:4443/api{}'.format(__api_version__),
    'test': 'http://api.salsa.skynet.unc.edu/{0[0]}.{0[1]}'.format(
        __api_version__),
    'production': 'https://api.skynet.unc.edu/{0[0]}.{0[1]}'.format(
        __api_version__),
}

warnings.filterwarnings('ignore', 'Unverified HTTPS request is being made')


class AttrDict(dict):
    """
    A dictionary-like object with obj.foo identical to obj['foo'] that also
    returns None for a missing key.
    """
    def __getattribute__(self, name):
        """
        Get dict keys as attributes

        :param name: attribute/key name

        :return: attribute/key value
        """
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            return dict.get(self, name)
    __getitem__ = __getattribute__

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __eq__(self, other):
        """
        Test equality of two dicts, compare floats to a reasonable precision

        :param other: value to compare with

        :return: True if objects are equal
        :rtype: bool
        """
        if not isinstance(other, dict):
            return False
        for key in set(tuple(self.keys()) + tuple(other.keys())):
            try:
                v1, v2 = self[key], other[key]
            except (AttributeError, KeyError):
                return False
            if isinstance(v1, float) and isinstance(v2, float):
                if abs(v1 - v2) > 1e-13:
                    return False
            elif v1 != v2:
                return False
        return True

    def __ne__(self, other):
        """
        Test inequality of two dicts, compare floats to a reasonable precision

        :param other: value to compare with

        :return: True if objects are not equal
        :rtype: bool
        """
        return not self.__eq__(other)

    def __copy__(self):
        """
        Return a new dict copy

        :rtype: AttrDict
        """
        return AttrDict(self)


def convert_json(value):
    """
    Convert dictionaries in a Python object to AttrDict; typical usage:
        res = convert_json(json.loads(json_string))

    :param value: Python object deserialized from a JSON string (a scalar,
        a list, or a dict)

    :return: input object with all dicts converted to AttrDict
    """
    if isinstance(value, dict):
        return AttrDict({key: convert_json(val) for key, val in value.items()})
    if isinstance(value, list):
        return [convert_json(item) for item in value]
    return value


def api_call(resource, method='get', params=None, server=None, api_key=None):
    """
    Call an API function

    :param str resource: API resource name, e.g. "obs", "users", or
        "groups/mygroup/users"
    :param str method: HTTP method: "get", "post", "put", or "delete"
    :param dict params: dictionary of request parameters
    :param str server: optional API server: "production", "test", "local", or an
        explicit URL; default is taken from ~/.skynet/api.cfg, assuming
        "production" if missing
    :param str api_key: optional API access token override; default: taken from
        ~/.skynet/API_KEY

    :return: result of the API call; raises the original API exception on error
    """
    if params is not None:
        # Convert request parameters to their appropriate string representation:
        # None -> empty string, date/time - ISO format, others -> as is
        params = {
            name: val.isoformat() if any(isinstance(val, t) for t in [
                datetime.datetime, datetime.date, datetime. time])
            else '' if val is None else val
            for name, val in params.items()}

    # Insert access token from ~/.skynet/API_KEY or its override into request
    # parameters
    if api_key is None:
        api_key = API_KEY
    if api_key:
        headers = {'Authentication-Token': api_key}

    else:
        headers = None

    # Make an API request to either the server from configuration or its
    # override
    kwargs = {}
    if params:
        if method.upper() in ('GET', 'HEAD', 'OPTIONS'):
            kwargs['params'] = params
        else:
            kwargs['data'] = params
    if headers:
        kwargs['headers'] = headers
    if server is None:
        server = DEFAULT_SERVER
    if server in BASE_URL:
        # if server != 'production':
        # Skip SSL certificate validation in debug mode; Skynet website
        # certificate temporarily does not provide the API access too
        kwargs['verify'] = False
        server = BASE_URL[server]

    r = requests.request(
        method.upper(), '{}/{}'.format(server, resource), **kwargs)

    if r.status_code != 200:
        # An error occurred
        exc = None
        if 'web2py_error' in r.headers:
            # Exception class passed in the response header, try to reconstruct
            exc_type = r.headers['web2py_error']
            # A hack: the original exception type object is restored from its
            # dotted package name (package.module.exctype) using the pickling
            # protocol 0 syntax; all importing of the necessary modules is
            # done by loads()
            try:
                mod, exc_type = exc_type.rsplit('.', 1)
            except ValueError:
                # No module name given? Assume builtin.
                mod = '__builtin__'
            try:
                exc = loads('c{}\n{}\np0\n.'.format(mod, exc_type).
                            encode('utf8'))
            except (AttributeError, ImportError):
                # Could not import from the given module, try skynet.api.errors
                pass

            if exc is None:
                try:
                    exc = getattr(errors, exc_type)
                except AttributeError:
                    # Not in skynet.api.errors as well, fall back to default
                    # exception type
                    pass

        if exc is None:
            # Could not deduce the exception class, raise the standard API error
            exc = errors.InternalApiError

        try:
            e = exc(r.text)
        except TypeError:
            # Exception instance requires other arguments than just the error
            # message
            raise errors.InternalApiError(r.text)
        raise e

    # On success, parse the result
    try:
        # Download-type request that returns a filename?
        try:
            fn = r.headers['Content-Disposition'].split(
                'attachment; filename=', 1)[1]
        except IndexError:
            fn = r.headers['Content-Disposition'].split(
                'inline; filename=', 1)[1]
    except (KeyError, IndexError):
        try:
            if r.headers['Content-Type'] == 'text/plain':
                # download/header-type request; return raw data as string
                return r.text
            # Normal request: deserialize from JSON and convert to AttrDict
            res = r.json()
            try:
                return convert_json(res)
            except (TypeError, ValueError):
                return res
        except Exception as e:
            raise errors.InternalApiError(str(e))
    else:
        # Result with an explicit file name and raw data; unquote filename and
        # return a pair (filename, data)
        try:
            return fn[fn.index('"') + 1:fn.rindex('"')].strip(), r.content
        except Exception as e:
            raise errors.InternalApiError(str(e))
