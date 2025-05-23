#!/usr/bin/env python

"""
Cancel observation(s)
"""

from __future__ import absolute_import, division, print_function

import argparse

from skynet.api.methods import obs as obs_api, radio_obs as radio_obs_api


def main():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        '-o', '--obs', metavar='ID,ID...',
        help='comma-separated list of observation IDs to cancel')
    parser.add_argument(
        '-u', '--user', metavar='ID/name,ID/name...',
        help='cancel observations submitted by specific user(s)')
    parser.add_argument(
        '-g', '--group', metavar='ID/name,ID/name...',
        help='cancel observations submitted by specific group(s)')
    parser.add_argument(
        '-s', '--scope', metavar='ID/name,ID/name...',
        help='cancel observations submitted to specific telescope/antenna(s)')
    parser.add_argument(
        '-n', '--antenna', metavar='ID/name,ID/name...',
        help='cancel observations submitted to specific radio telescope(s)')
    parser.add_argument(
        '-a', '--after', metavar='yyyy-mm-ddThh:mm:ss.s',
        help='cancel observations submitted after the given UTC')
    parser.add_argument(
        '-r', '--radio', action='store_true', help='cancel radio observations')
    parser.add_argument(
        '-b', '--before', metavar='yyyy-mm-ddThh:mm:ss.s',
        help='cancel observations submitted before the given UTC')

    args = parser.parse_args()
    kw = {key: val for key, val in getattr(args, '_get_kwargs')()
          if val is not None}  # suppress protected member warn
    obs_ids = kw.pop('obs', '')
    optical = not kw.pop('radio', False)

    if optical:
        if obs_ids:
            obs_ids = obs_ids.split(',')
        else:
            # Get optical observation IDs:
            #     skynet.api.methods.obs.query(
            #         scope='Prompt1,Prompt3',
            #         after=datetime.datetime.utcnow() -
            #         datetime.timedelta(days=1))
            # Supported keywords: obs, user, group, scope, after, before
            obs_ids = obs_api.query(**kw)

        # Cancel observations one by one by setting state to "canceled"
        # noinspection PyTypeChecker
        for obs in obs_ids:
            obs_api.update(obs, state='canceled')
    else:
        if obs_ids:
            obs_ids = obs_ids.split(',')
        else:
            # Get radio observation IDs:
            #     skynet.api.methods.radio_obs.query(
            #         antenna='GreenBank-20',
            #         after=datetime.datetime.utcnow() -
            #         datetime.timedelta(days=1))
            # Supported keywords: obs, user, group, antenna, after, before
            obs_ids = radio_obs_api.query(**kw)

        # Cancel observations one by one by setting state to "canceled"
        # noinspection PyTypeChecker
        for obs in obs_ids:
            radio_obs_api.update(obs, state='canceled')


if __name__ == '__main__':
    main()
