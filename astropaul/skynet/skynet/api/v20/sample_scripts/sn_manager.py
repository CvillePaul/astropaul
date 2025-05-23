#!/usr/bin/env python

"""
Manage supernovae search observations
"""

from __future__ import absolute_import, division, print_function

import os.path
import sys
import time
import argparse
import json
import datetime
from math import acos, asin, atan, cos, floor, sin, pi, tan

from skynet.api.client import convert_json
from skynet.api.methods import (
    exps as exp_api, download as download_api, obs as obs_api,
    scopes as scope_api)
from skynet.api.errors import ResourceNotFoundError


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""\
The list of targets is stored in a JSON-formatted file as follows:

[
    {
        "name": "unique name of the first target",
        "raHours": RA of the target in decimal hours or as "HH:MM:SS.SS",
        "decDegs": Dec of the target in decimal degrees or as "+DD:MM:SS.S",
        "telescopes": "name1, ID2..." or ["name1", ID2] - optional,
        "timeAccountId": time account ID - optional,
        "objectName": "target name" - optional,
        "maxSun": maximum Sun elevation in decimal degrees - optional,
        "priority": relative priority of observation (integer >= 0) - optional,
        ... a number of other observation fields (see API docs)
        "exps": [
            {
                "expLength": exposure length in seconds,
                "filterRequested": "filter name" or integer ID,
                "telescopeRequested": "name" or ID - optional
            },
            {
                ... next exposure spec
            },
            ...
        ]
    },
    {
        ... next target spec
    },
    ...
]""")
    parser.add_argument(
        '-a', '--time-account', metavar='ID', type=int, help='time account ID')
    parser.add_argument(
        '-d', '--delay', type=float, default=60, metavar='seconds',
        help='how long to delay between loops, in seconds (if --repeat has '
             'been specified)')
    parser.add_argument(
        '-i', '--force-int', action='store_true',
        help='always download images in the 16-bit unsigned format')
    parser.add_argument(
        '-r', '--repeat', action='store_true',
        help='repeat script continuously (default: run once)')
    parser.add_argument(
        '-t', '--telescopes', metavar='scope1,scope2...',
        help='comma-separated list of telescope IDs or names, first is primary')
    parser.add_argument(
        '-u', '--reduce', action='store_true', help='reduce images on download')
    parser.add_argument('target_list', help='name of target list file')

    args = parser.parse_args()

    if args.repeat:
        print('Running continuously, press Ctrl-C to terminate')

    while True:
        try:
            run(args)
        except (KeyboardInterrupt, SystemExit):
            break
        except Exception as e:
            print(', '.join(str(arg) for arg in e.args) if e.args else e,
                  file=sys.stderr)
        if not args.repeat:
            break

        try:
            time.sleep(args.delay)
        except (KeyboardInterrupt, SystemExit):
            break


target_file_timestamp = None

# Load the list of already downloaded exposures
downloaded_exps_filename = os.path.expanduser('~/.skynet/.sn_manager.exps')
downloaded_exps = []
try:
    with open(downloaded_exps_filename, 'rU') as _f:
        for _exp_id in _f.read().splitlines():
            try:
                downloaded_exps.append(int(_exp_id))
            except ValueError:
                pass
except IOError:
    pass

failed_exps = []


def id_and_name(exp, id_attr, name_attr):
    """
    Helper function for identical_exps() that returns ID and name (both may be
    None) of the given field, e.g. exp.filterIdRequested and exp.filterId;
    works for all forms, including
        exp.filterIdRequested = id,
        exp.filterRequested = id,
        exp.filterRequested = 'id',
        exp.filterRequested = 'bar',
        exp.filterRequested = {'id':id}, {'name':'bar'}, {'id':id,'name':'bar'}

    :param skynet.api.client.AttrDict exp: exposure spec
    :param str id_attr: name of ID attribute (e.g. 'filterIdRequested')
    :param str name_attr: name attribute (e.g. 'filterRequested')

    :return: (id, name)
    :rtype: tuple[int | None, str | None]
    """
    i, name = getattr(exp, id_attr), getattr(exp, name_attr)
    try:
        if name.id is not None:
            i = name.id
        name = name.name
    except AttributeError:
        try:
            i = int(name)
            name = None
        except (TypeError, ValueError):
            pass
    return i, name


def identical_field(e1, e2, id_attr, name_attr):
    """
    Return True if two exposure have the same value of the given parameter;
    used by identical_exps()

    :param skynet.api.client.AttrDict e1: first exposure spec
    :param skynet.api.client.AttrDict e2: second exposure spec:
    :param str id_attr: name of ID attribute (e.g. 'filterIdRequested')
    :param str name_attr: name attribute (e.g. 'filterRequested')

    :rtype: bool
    """
    id1, name1 = id_and_name(e1, id_attr, name_attr)
    id2, name2 = id_and_name(e2, id_attr, name_attr)
    return (id1 is not None and id1 == id2 or
            name1 is not None and name1 == name2)


def identical_exps(e1, e2):
    """
    Return True if two exposure have the same length, filter, and telescope;
    used to check whether a new exposure should be submitted for an existing obs

    :param skynet.api.client.AttrDict e1: first exposure spec
    :param skynet.api.client.AttrDict e2: second exposure spec:

    :rtype: bool
    """
    return abs(e1.expLength - e2.expLength) < 1e-6 and \
        identical_field(e1, e2, 'filterIdRequested', 'filterRequested') and \
        identical_field(e1, e2, 'teleIdRequested', 'telescopeRequested')


def get_rise_set_time(sunrise, d, lat, lon, h=-12):
    """
    Calculate the sunrise/sunset time at the given site

    :param bool sunrise: True for sunrise, False for sunset
    :param datetime.date | datetime.datetime d: UTC date
    :param float lat: latitude in degrees (+N)
    :param float lon: longitude in degrees (+E)
    :param float h: optional Sun elevation in degrees

    :return: UTC time of sunrise/sunset in hours or None if there's no
        sunrise/sunset on the given date
    :rtype: datetime.datetime | None
    """
    if isinstance(d, datetime.datetime):
        d = d.date()
    n = int((d - datetime.date(d.year, 1, 1)).total_seconds()//86400 + 1.5)
    lon_h = lon/15
    t = n + (18 - 12*bool(sunrise) - lon_h)/24
    m = 0.9856*t - 3.289
    l = (m + 1.916*sin(m*pi/180) + 0.020*sin(2*m*pi/180) + 282.634) % 360
    ra = (atan(0.91764*tan(l*pi/180))*180/pi) % 360
    ra += floor(l/90)*90 - floor(ra/90)*90
    ra /= 15
    sin_dec = 0.39782*sin(l*pi/180)
    cos_dec = cos(asin(sin_dec))
    cos_h = (cos((90 - h)*pi/180) - sin_dec*sin(lat*pi/180))/(
        cos_dec*cos(lat*pi/180))
    if abs(cos_h) > 1:
        return None
    h = acos(cos_h)*12/pi
    return datetime.datetime.combine(d, datetime.time()) + datetime.timedelta(
        hours=(h*(1 - 2*sunrise) + ra - (0.06571*t) - 6.622 - lon_h) % 24)


def run(args):
    """
    Do a single SN observation management run

    :param argparse.Namespace args: scripr arguments

    :rtype: None
    """
    global target_file_timestamp

    # Load the list of targets if changed
    ts = os.stat(args.target_list).st_mtime
    if target_file_timestamp == ts:
        return

    with open(args.target_list, 'rU') as f:
        try:
            target_list = convert_json(json.load(f))
        except Exception as e:
            target_list = None
            raise Exception('{}: {}'.format(
                os.path.basename(args.target_list),
                ', '.join(str(arg) for arg in e.args) if e.args else e))

    # Manage each observation in the target list
    for target in target_list:
        name = target.name
        exps = target.pop('exps', [])
        if not exps:
            print('WARNING. Empty exposure list for "{}" in "{}"'.format(
                name, args.target_list), file=sys.stderr)

        if args.time_account is not None and 'timeAccountId' not in target:
            target.timeAccountId = args.time_account
        if args.telescopes is not None and 'telescopes' not in target:
            target.telescopes = args.telescopes

        try:
            obs = obs_api.get(name)
        except ResourceNotFoundError:
            # Observation for the target does not exist, submit an empty one
            print('Submitting a new observation for "{}"'.format(name))
            obs = obs_api.add(**target)
            print('Obs ID: {:d}, telescopes: {}'.format(
                obs.id, ', '.join(scope.name for scope in obs.telescopes)))
        else:
            del obs.exps

            # Observation parameters could have changed, update as needed and
            # refresh
            old_obs = obs
            obs = obs_api.update(obs.id, **target)
            del obs.exps
            if obs != old_obs:
                print('Updated observation for "{}" ({:d})'.format(
                    name, obs.id))
                for key in sorted(set(
                        tuple(old_obs.keys()) + tuple(obs.keys()))):
                    if key in old_obs:
                        if key in obs:
                            v1, v2 = old_obs[key], obs[key]
                            if isinstance(v1, float) and isinstance(v2, float):
                                if abs(v1 - v2) > 1e-13:
                                    print('{}: {} -> {}, diff: {}'.format(
                                        key, v1, v2, v2 - v1))
                            elif v1 != v2:
                                print('{}: {} -> {}'.format(key, v1, v2))
                        else:
                            print('{}: {} -> MISSING'.format(key, old_obs[key]))
                    else:
                        print('{}: MISSING ->'.format(key, obs[key]))

        # (Re-)submit missing exposures for the current observation
        submitted = []
        for exp_spec in exps:
            # Count requested and active exposures with the same parameters
            n_requested = len([
                exp for exp in exps if identical_exps(exp, exp_spec)])
            if n_requested <= len([
                    1 for exp in exp_api.query(
                        obs=obs.id, state='ready',
                        include='id,expLength,filterIdRequested,'
                        'filterRequested,teleIdRequested,telescopeRequested')
                    if identical_exps(exp, exp_spec)]):
                # Enough active exposures, do nothing
                continue

            # Not enough active exposures with the given parameters; submit
            # during daytime or if no such exposures were successfully taken
            # since the last sunset
            scope_id, scope_name = id_and_name(
                exp_spec, 'teleIdRequested', 'telescopeRequested')
            if scope_id or scope_name:
                # Explicit telescope assignment, get from exp
                if scope_id:
                    try:
                        site = scope_api.get(scope_id).site
                    except IndexError:
                        print('WARNING. Unknown telescope ID: {}'.format(
                            scope_id), file=sys.stderr)
                        continue
                else:
                    try:
                        site = scope_api.get(scope_name).site
                    except IndexError:
                        print('WARNING. Unknown telescope: "{}"'.format(
                            scope_name), file=sys.stderr)
                        continue
            else:
                # No explicit assignment, get from obs; choose any site if the
                # obs is requested for multiple sites
                try:
                    site = list(set([scope.site
                                     for scope in obs.telescopes]))[0]
                except IndexError:
                    # No scopes set yet, this means that it is the first
                    # submission, and it should proceed
                    site = None
            if site:
                t = datetime.datetime.utcnow()
                t_rise = get_rise_set_time(True, t, site.latDegs, site.lngDegs)
                if t_rise is None:
                    t_rise = t
                t_set = get_rise_set_time(False, t, site.latDegs, site.lngDegs)
                if t_set is None:
                    t_set = t
                if (t_rise > t_set and t_set <= t <= t_rise or
                        t_rise <= t_set and (
                            t_set <= t <= t_rise + datetime.timedelta(days=1) or
                            t_set - datetime.timedelta(days=1) <= t <= t_rise))\
                        and len([
                            1 for exp in exp_api.query(
                                obs=obs.id, state='completed,archived',
                                after=t_set - datetime.timedelta(
                                    days=1 if t <= t_set else 0),
                                include='id,expLength,filterIdRequested,'
                                'filterRequested,teleIdRequested,'
                                'telescopeRequested')
                            if identical_exps(exp, exp_spec)]
                        ) >= n_requested:
                    # Night time, observations still going on, and the same
                    # exposure was completed during this night
                    continue
            exp = exp_api.add(obs.id, **exp_spec)
            submitted.append(str(exp.id))
        if submitted:
            print('Submitted {:d} exposure(s) for "{}" (obs {:d}): '
                  '{}'.format(
                      len(submitted), obs.name, obs.id,
                      ', '.join(submitted)))

        # Download all exposures not yet downloaded
        for exp in exp_api.query(
                obs=obs.id, state='archived', wcs_present=1):
            if exp.id in downloaded_exps:
                continue

            try:
                filename, data = download_api.fits(
                    image='r{:d}'.format(exp.id), reduce=int(args.reduce),
                    force_int=int(args.force_int))
            except Exception as e:
                if exp.id not in failed_exps:
                    failed_exps.append(exp.id)
                    print('WARNING. Error downloading exposure ID {:d}: '
                          '{}'.format(
                              exp.id,
                              ', '.join(str(arg) for arg in e.args)
                              if e.args else e), file=sys.stderr)
            else:
                print('Downloaded exposure ID {:d} ({})'.format(
                    exp.id, filename))
                try:
                    with open(filename, 'wb') as f:
                        f.write(data)
                except Exception as e:
                    if exp.id not in failed_exps:
                        failed_exps.append(exp.id)
                        print('WARNING. Error saving image "{}": {}'.format(
                            filename,
                            ', '.join(str(arg) for arg in e.args) if e.args
                            else e), file=sys.stderr)
                else:
                    downloaded_exps.append(exp.id)
                    try:
                        failed_exps.remove(exp.id)
                    except ValueError:
                        pass

                    try:
                        with open(downloaded_exps_filename, 'at') as f:
                            f.write('{:d}\n'.format(exp.id))
                    except Exception as e:
                        print('WARNING. Error updating the list of downloaded '
                              'images for exposure ID {:d}: {}\n'
                              'Expect repeated downloads of the same file on '
                              'future runs.'.format(
                                  exp.id,
                                  ', '.join(str(arg) for arg in e.args)
                                  if e.args else e), file=sys.stderr)

    # Target updated successfully
    target_file_timestamp = ts


if __name__ == '__main__':
    main()
