#!/usr/bin/env python
"""
SkyNet API test script

Usage: test_api.py method resource [param=value ...]
"""

from __future__ import absolute_import, division, print_function
import sys
import pprint
from skynet.api.client import api_call


def main():
    method, resource = sys.argv[1:3]

    params = {}
    if resource == 'obs' and method.lower() == 'post':
        params.update(dict(name='test', raHours=0, decDegs=0))

    # noinspection PyBroadException
    server = username = None
    for arg in sys.argv[3:]:
        name, val = arg.split('=', 1)
        if name == 'SERVER':
            server = val
        elif name == 'USER':
            username = val
        else:
            params[name] = val

    if username:
        # Specific user, get token from db (server side only)
        try:
            # noinspection PyUnresolvedReferences
            from skynet.db.model import DatabaseSession, OauthCode, User
        except ImportError:
            raise RuntimeError(
                'USER=... is allowed on the SkyNet server side only')
        adb = DatabaseSession()
        user = adb.query(User)
        try:
            user = user.get(int(username))
        except ValueError:
            user = user.filter(User.username == username).one_or_none()
        if user is None:
            raise RuntimeError('Unknown user: "{}"'.format(username))
        user_id, username = user.id, user.username
        try:
            api_key = adb.query(OauthCode).filter(
                OauthCode.userId == user_id).first().tokens[0].accessToken
        except AttributeError:
            raise RuntimeError('User "{}" ({:d}) has no API key'.format(
                username, user_id))
    else:
        # Default user, assume token is in ~/.skynet
        api_key = None

    res = api_call(resource, method, params, server, api_key)
    if isinstance(res, tuple):
        # (filename, data) tuple for a download request
        if res[0]:
            filename = str(res[0])
            print('Downloading {} ({} byte(s))'.format(
                filename, len(res[1])))
            with open(filename, 'wb') as f:
                f.write(res[1])
        else:
            print(res[1])
    elif isinstance(res, str) or isinstance(res, type(u'')) or \
            isinstance(res, list):
        print(res)
    else:
        pprint.PrettyPrinter(indent=4).pprint(res)


if __name__ == '__main__':
    main()
