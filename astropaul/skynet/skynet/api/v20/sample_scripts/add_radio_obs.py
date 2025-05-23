#!/usr/bin/env python

"""
Submit a radio observation
"""

from __future__ import absolute_import, division, print_function

import sys
import argparse

from skynet.api.methods import radio_obs as radio_obs_api


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog='Prints obs ID and credits charged on successful submission')

    parser.add_argument(
        'type', metavar='track|map|daisy|onoff', help='observation type')
    parser.add_argument('name', metavar='"name"', help='observation name')
    parser.add_argument(
        'duration', type=float, metavar='duration', help='duration in seconds')
    parser.add_argument(
        'mode', metavar='lowres|highres|pulsar', help='observation mode')
    parser.add_argument(
        '--ra', metavar='degrees|DD:MM:SS.S',
        help='RA of the target object; requires --dec and is mutually '
        'exclusive with --lat, --lon, --az, and --el')
    parser.add_argument(
        '--dec', metavar='degrees|+DD:MM:SS.S',
        help='Dec of the target object; requires --ra and is mutually '
        'exclusive with --lat, --lon, --az, and --el')
    parser.add_argument(
        '--lon', metavar='degrees|DD:MM:SS.S',
        help='Galactic longitude of the target object; requires --lat and is '
        'mutually exclusive with --ra, --dec, --az, and --el')
    parser.add_argument(
        '--lat', metavar='degrees|+DD:MM:SS.S',
        help='Galactic latitude of the target object; requires --lon and is '
        'mutually exclusive with --ra, --dec, --az, and --el')
    parser.add_argument(
        '--az', metavar='degrees|+DD:MM:SS.S',
        help='Azimuth of the target object; requires --el and is mutually '
        'exclusive with --ra, --dec, --lat, and --lon')
    parser.add_argument(
        '--el', metavar='degrees|+DD:MM:SS.S',
        help='Elevation of the target object; requires --az and is mutually '
        'exclusive with --ra, --dec, --lat, and --lon')
    parser.add_argument(
        '--object-type', metavar='sidereal|planet|asteroid|comet',
        dest='centerObjectType', help='target object type')
    parser.add_argument(
        '--object-name', metavar='"name"', dest='centerObjectName',
        help='target object name if --obj-type is not sidereal')
    parser.add_argument(
        '--start-after', metavar='yyyy-mm-ddThh:mm:ss.s', dest='startAfter',
        help='start observation after the given UTC')
    parser.add_argument(
        '--end-before', metavar='yyyy-mm-ddThh:mm:ss.s', dest='endBefore',
        help='do NOT start observation after the given UTC')
    parser.add_argument(
        '--priority', metavar='N', type=int, help='priority of observation')
    parser.add_argument(
        '--channels', metavar='N', type=int, dest='nChannel',
        help='number of channels')
    parser.add_argument(
        '--center-freq', metavar='MHz', type=float, dest='centerFrequency',
        help='central frequency')
    parser.add_argument(
        '--secondary-freq', metavar='MHz', type=float,
        dest='secondaryFrequency', help='secondary frequency')
    parser.add_argument(
        '--integration', metavar='s', type=float, dest='integrationTime',
        help='integration time')
    parser.add_argument(
        '--filter', metavar='FILTER', dest='radioFilter',
        help='radio filter name or ID')
    parser.add_argument(
        '--min-el', metavar='degrees|+DD:MM:SS.S', dest='minEl',
        help='minimum allowed elevation of the target object')
    parser.add_argument(
        '--min-solar-sep', metavar='degrees|DD:MM:SS', dest='minSolarSep',
        help='minimum allowed solar separation of the target object')
    parser.add_argument(
        '-a', '--time-account', metavar='ID', type=int, dest='timeAccountId',
        help='time account ID')
    parser.add_argument(
        '-t', '--antennas', metavar='antenna1,antenna2...',
        help='comma-separated list of antenna IDs or names')
    parser.add_argument(
        '--too', metavar='"justification"',
        help='request TOO with the given justification')
    parser.add_argument(
        '-r', '--repeat', metavar='N', type=int,
        help='type=track|onoff: number of repeats')
    parser.add_argument(
        '--start-hoffset', metavar='degrees|+DD:MM:SS.S', dest='fixedOffsetH',
        help='type=track: starting RA/Az/Lng offset')
    parser.add_argument(
        '--start-voffset', metavar='degrees|+DD:MM:SS.S', dest='fixedOffsetV',
        help='type=track: starting Dec/El/Lat offset')
    parser.add_argument(
        '--end-hoffset', metavar='degrees|+DD:MM:SS.S', dest='endOffsetH',
        help='type=track: ending RA/Az/Lng offset')
    parser.add_argument(
        '--end-voffset', metavar='degrees|+DD:MM:SS.S', dest='endOffsetV',
        help='type=track: ending Dec/El/Lat offset')
    parser.add_argument(
        '--map-type', metavar='declat|ralong', dest='mapType',
        help='type=map: map direction (required)')
    parser.add_argument(
        '--hlength', metavar='degrees|DD:MM:SS.S', dest='hLength',
        help='type=map: RA/Az/Lng map size (required)')
    parser.add_argument(
        '--vlength', metavar='degrees|DD:MM:SS.S', dest='vLength',
        help='type=map: Dec/El/Lat map size (required)')
    parser.add_argument(
        '--delta', metavar='degrees|DD:MM:SS.S',
        help='type=map: gap between sweeps (required)')
    parser.add_argument(
        '--unidirectional', metavar='0|1', type=int,
        help='type=map: unidirectional map')
    parser.add_argument(
        '--radius', type=float, metavar='arcminutes',
        help='type=daisy: daisy radius (required)')
    parser.add_argument(
        '--period', type=float, metavar='s', dest='radialPeriod',
        help='type=daisy: radial period (required)')
    parser.add_argument(
        '--radial-phase', metavar='degrees|+DD:MM:SS.S', dest='radialPhase',
        help='type=daisy: radial phase')
    parser.add_argument(
        '--rotation-phase', metavar='degrees|+DD:MM:SS.S', dest='rotationPhase',
        help='type=daisy: rotation phase')
    parser.add_argument(
        '--hoffset', metavar='degrees|+DD:MM:SS.S', dest='referenceOffsetH',
        help='type=onoff: RA/Az/Lng reference offset')
    parser.add_argument(
        '--voffset', metavar='degrees|+DD:MM:SS.S', dest='referenceOffsetV',
        help='type=onoff: Dec/El/Lat reference offset')
    parser.add_argument(
        '--reverse', metavar='0|1', type=int,
        help='type=onoff: reverse direction')

    args = parser.parse_args()
    kw = {key: val for key, val in getattr(args, '_get_kwargs')()
          if val is not None}  # suppress protected member warning

    # Most arguments are named exactly as skynet.api.methods.radio_obs.add()
    # expects; handle a few special cases
    kw['obsType'] = args.type
    kw['recordType'] = args.mode

    if args.ra or args.dec:
        if not (args.ra and args.dec):
            print('Both --ra and --dec must be given', file=sys.stderr)
            sys.exit(1)
        if args.lat or args.lon or args.az or args.el:
            print('Equatorial coordinates are mutually exclusive with galactic '
                  'and horizontal coordinates', file=sys.stderr)
            sys.exit(1)
        kw['centerH'], kw['centerV'] = args.ra, args.dec
        kw['coordType'] = 'RA_DEC_COORD'
    elif args.lon or args.lat:
        if not (args.lon and args.lat):
            print('Both --lon and --lat must be given', file=sys.stderr)
            sys.exit(1)
        if args.az or args.el:
            print('Galactic coordinates are mutually exclusive with equatorial '
                  'and horizontal coordinates', file=sys.stderr)
            sys.exit(1)
        kw['centerH'], kw['centerV'] = args.lon, args.lat
        kw['coordType'] = 'LNG_LAT_COORD'
    elif args.az or args.el:
        if not (args.az and args.el):
            print('Both --az and --el must be given', file=sys.stderr)
            sys.exit(1)
        kw['centerH'], kw['centerV'] = args.az, args.el
        kw['coordType'] = 'AZ_EL_COORD'
    else:
        print('Missing target object coordinates', file=sys.stderr)
        sys.exit(1)

    if args.too:
        kw['isToo'] = True
        kw['tooJustification'] = args.too

    for key in ('type', 'mode', 'ra', 'dec', 'lat', 'lon', 'az', 'el', 'too'):
        try:
            del kw[key]
        except KeyError:
            pass

    # Submit observation and print its ID and credits charged
    obs = radio_obs_api.add(**kw)
    print(obs.id, obs.creditsCharged)


if __name__ == '__main__':
    main()
