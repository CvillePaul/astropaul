#!/usr/bin/env python

"""
Submit an observation
"""

from __future__ import absolute_import, division, print_function

import sys
import argparse
import json

from skynet.api.methods import obs as obs_api


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog='Prints serialized Observation object on successful submission')

    parser.add_argument('name', metavar='"name"', help='observation name')
    parser.add_argument(
        '--ra', metavar='hours|HH:MM:SS.S', dest='raHours',
        help='RA of the target object')
    parser.add_argument(
        '--dec', metavar='degrees|+DD:MM:SS.S', dest='decDegs',
        help='Dec of the target object')
    parser.add_argument(
        '--exp-lengths', metavar='t1,t2,...,tn',
        help='comma-separated list of exposure lengths')
    parser.add_argument(
        '--exp-filters', metavar='f1,f2,...,fn',
        help='comma-separated list of filter IDs/names')
    parser.add_argument(
        '--exp-scopes', metavar='s1,s2,...,sn',
        help='comma-separated list of telescope IDs/names for each exposure; '
        'leave blank to use defaults from --telescopes, e.g. '
        'Prompt1,12,,Propmt3,')
    parser.add_argument(
        '--priority', metavar='N', type=int, help='priority of observation')
    parser.add_argument(
        '-a', '--time-account', metavar='ID', type=int, dest='timeAccountId',
        help='time account ID')
    parser.add_argument(
        '-t', '--telescopes', metavar='scope1,scope2...',
        help='comma-separated list of telescope IDs or names')
    parser.add_argument(
        '-e', '--efficiency', metavar='scope',
        help='assume exp lengths for this efficiency')
    parser.add_argument(
        '--too', metavar='"justification"',
        help='request TOO with the given justification')
    parser.add_argument(
        '--min-el', type=float, metavar='degrees', dest='minEl',
        help='minimum elevation')
    parser.add_argument(
        '--max-sun', type=float, metavar='degrees', dest='maxSun',
        help='maximum Sun elevation')
    parser.add_argument(
        '--object-name', metavar='"name"', dest='objectName',
        help='target object name')
    parser.add_argument(
        '--object-type', metavar='sidereal|planet|asteroid|comet|nonsidereal',
        dest='objectType', help='target object type')
    parser.add_argument(
        '--object-dist', type=float, metavar='AU', dest='objectDist',
        help='distance to solar system object, in AU')
    parser.add_argument(
        '--cancel-after', metavar='yyyy-mm-ddThh:mm:ss.s',
        dest='cancelAfterUtc',
        help='automatically cancel observation after the given UTC')
    parser.add_argument(
        '--next-exp-after', metavar='yyyy-mm-ddThh:mm:ss.s',
        dest='nextExpStartAfterUtc',
        help='start next exposure after the given UTC')
    parser.add_argument(
        '--dither', action='store_true', dest='ditherEnabled',
        help='enable dithering')
    parser.add_argument(
        '--dither-x', type=int, metavar='pixels', dest='ditherXSize',
        help='horizontal dithering size')
    parser.add_argument(
        '--dither-y', type=int, metavar='pixels', dest='ditherYSize',
        help='vertical dithering size')
    parser.add_argument(
        '--dither-spacing', type=float, metavar='arcsecs',
        dest='ditherSpacingArcsecs', help='dither spacing')
    parser.add_argument(
        '--tracking', metavar='track_target|lock_field', dest='targetTracking',
        help='target tracking mode')
    parser.add_argument(
        '--field-lock-at', metavar='yyyy-mm-ddThh:mm:ss.s',
        dest='fieldLockUtc',
        help='lock field at the given UTC whith "--tracking lock_field"')
    parser.add_argument(
        '--trigger-repoint', type=float, metavar='arcmins',
        help='repoint threshold')
    parser.add_argument(
        '--point-ahead', type=float, metavar='secs', help='enable point-ahead')
    parser.add_argument(
        '--ra-offset', type=float, metavar='arcmins',
        dest='constantRaOffsetArcmins', help='constant RA offset')
    parser.add_argument(
        '--dec-offset', type=float, metavar='arcmins',
        dest='constantDecOffsetArcmins', help='constant Dec offset')
    parser.add_argument(
        '--moon-sep', type=float, metavar='degrees', dest='minMoonSepDegs',
        help='minimum Moon separation')
    parser.add_argument(
        '--rbi-limit', type=float, dest='rbiFractionAvgBkgLimit',
        help='average background level fraction to trigger RBI protection')

    # Most arguments are named exactly as skynet.api.methods.obs.add() expects;
    # handle a few special cases
    args = parser.parse_args()
    kw = {key: val for key, val in getattr(args, '_get_kwargs')()
          if val is not None}  # suppress protected member warn
    if args.too:
        kw['isToo'] = True
        kw['tooJustification'] = args.too
    if args.trigger_repoint:
        kw['triggerRepointEnabled'] = True
        kw['triggerRepointArcmins'] = args.trigger_repoint
    if args.point_ahead:
        kw['pointAheadEnabled'] = True
        kw['pointAheadSecs'] = args.point_ahead
    for key in ('too', 'trigger_repoint', 'point_ahead', 'exp_lengths',
                'exp_filters', 'exp_scopes'):
        try:
            del kw[key]
        except KeyError:
            pass

    # Extract exposure specs
    if args.exp_lengths:
        exp_lengths = [float(item) for item in args.exp_lengths.split(',')]
    else:
        exp_lengths = []
    if args.exp_filters:
        exp_filters = [item.strip() for item in args.exp_filters.split(',')]
    else:
        exp_filters = []
    if args.exp_scopes:
        exp_scopes = [item.strip() for item in args.exp_scopes.split(',')]
    else:
        exp_scopes = []
    n_exp = max(len(exp_lengths), len(exp_filters), len(exp_scopes))
    if n_exp:
        if len(exp_lengths) < n_exp:
            if not exp_lengths:
                print('Missing --exp-lengths', file=sys.stderr)
                sys.exit(1)
            exp_lengths += [exp_lengths[-1]]*(n_exp - len(exp_lengths))
        if len(exp_filters) < n_exp:
            if not exp_filters:
                print('Missing --exp-filters', file=sys.stderr)
                sys.exit(1)
            exp_filters += [exp_filters[-1]]*(n_exp - len(exp_filters))
        if len(exp_scopes) < n_exp:
            if not exp_scopes:
                exp_scopes = [None]*n_exp
            else:
                exp_scopes += [exp_scopes[-1]]*(n_exp - len(exp_scopes))

    # Submit observation:
    #     obs = obs_api.add(name='name', raHours='1:2:3', decDegs=-1.23, ...)
    # Supported keywords: name, raHours, decDegs, priority, timeAccountId,
    # telescopes, efficiency, isToo, tooJustification, minEl, maxSun,
    # objectName, objectType, objectDist, cancelAfterUtc, nextExpStartAfterUtc,
    # ditherEnabled, ditherXSize, ditherYSize, ditherSpacingArcsecs,
    # targetTracking, fieldLockUtc, triggerRepointEnabled,
    # triggerRepointArcmins, pointAheadEnabled, pointAheadSecs,
    # constantRaOffsetArcmins, constantDecOffsetArcmins, minMoonSepDegs,
    # rbiFractionAvgBkgLimit
    obs = obs_api.add(
        exps=json.dumps([
            dict(expLength=l, filterRequested=f, telescopeRequested=t)
            for l, f, t in zip(exp_lengths, exp_filters, exp_scopes)
        ]),
        **kw)
    print(obs)


if __name__ == '__main__':
    main()
