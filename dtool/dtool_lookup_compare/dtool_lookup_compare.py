#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
#
# dtool_lookup_compare.py
#
# Copyright (C) 2021 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Compare local repository of dtool datasets against registered datasets on lookup server.
'''
import logging
import argcomplete
import argparse

import json
import dtool_lookup_api as dl

from datetime import datetime

import dtoolcore
from dtool_cli.cli import CONFIG_PATH
from dtool_info.utils import sizeof_fmt, date_fmt


def _make_marker(d):
    """Mark everything for comparison."""
    if isinstance(d, list):
        return [_make_marker(e) for e in d]
    elif isinstance(d, dict):
        return {k: _make_marker(v) for k, v in d.items()}
    else:
        return True

def _compare(source, target, marker):
    """Compare source and target partially, as marked by marker."""
    logger = logging.getLogger(__name__)
    if isinstance(marker, dict):
        for k, v in marker.items():
            if k not in source:
                logger.error("{} not in source '{}'.".format(k, source))
                return False
            if k not in target:
                logger.error("{} not in target '{}'.".format(k, source))
                return False

            logger.debug("Descending into sub-tree '{}' of '{}'.".format(
                source[k], source))
            # descend
            if not _compare(source[k], target[k], v):
                return False  # one failed comparison suffices

    elif isinstance(marker, list):  # source, target and marker must have same length
        logger.debug("Branching into element wise sub-trees of '{}'.".format(
            source))
        for s, t, m in zip(source, target, marker):
            if not _compare(s, t, m):
                return False  # one failed comparison suffices
    else:  # arrived at leaf, comparison desired?
        if marker is not False:  # yes
            logger.debug("Comparing '{}' == '{}' -> {}.".format(
                source, target, source == target))
            return source == target

    # comparison either not desired or successfull for all elements
    return True


def _compare_nested(source, target, marker=None):
    """Compare source and target partially, as marked by marker. If marker is None, then compare everything."""
    if not marker:
        marker = _make_marker(source)
    return _compare(source, target, marker)


def _forward_compare(source, target, marker=None):
    """One-way-compare two dicts of nested dict and categorize into 'equal', 'differing' and 'missing'."""
    missing = dict()
    differing = dict()
    equal = dict()

    for k, sd in source.items():
        if k in target:
            td = target[k]
            is_equal = _compare_nested(sd, td, marker)
            if is_equal:
                equal[k] = (sd, td)
            else:
                differing[k] = (sd, td)
        else:
            missing[k] = sd

    return equal, differing, missing


def _ds_list_to_dict(l):
    """Convert list of dataset metadata entries to dict with UUIDs as keys."""
    return {e['uuid']: e for e in l}


def _direct_list(base_uri, config_path=CONFIG_PATH):
    """Directly list all datasets at base_uri via suitable storage broker."""
    base_uri = dtoolcore.utils.sanitise_uri(base_uri)
    StorageBroker = dtoolcore._get_storage_broker(base_uri, config_path)
    info = []

    for uri in StorageBroker.list_dataset_uris(base_uri, config_path):
        admin_metadata = dtoolcore._admin_metadata_from_uri(uri, config_path)
        info.append(admin_metadata)

    by_name = sorted(info, key=lambda d: d['name'])
    return by_name


def _lookup_list(query={}):
    """List all datasets registered at lookup servr, filtered by query."""
    res = dl.query(query)
    by_name = sorted(res, key=lambda d: d['name'])
    return by_name


def _print_nested(l):
    """Print a nested structure, i.e list or dict."""
    print(json.dumps(l, indent=4))


def _print_dataset_list(l, verbose=True):
    """Print a list of dataset metadata entries."""
    if not verbose:
        l = [e[0]['uuid'] if isinstance(e, tuple) else e['uuid'] for e in l]
    _print_nested(l)


def _compare_dataset_lists(source, target, marker=None):
    """One-way-compare two lists of dataset metadata entries."""
    s = _ds_list_to_dict(source)
    t = _ds_list_to_dict(target)
    equal, differing, missing = _forward_compare(s, t, marker)
    return list(equal.values()), list(differing.values()), list(missing.values())


def compare_dataset_lists(source, target,
                          marker={'uuid': True, 'name': True, 'type': True, 'created_at': True, 'frozen_at': True},
                          verbose=True, terse=False,
                          print_equal=True,
                          print_differing=True,
                          print_missing=True):
    """Compare source and target dataset metadata lists by fields set True within marker."""
    equal, differing, missing = _compare_dataset_lists(source, target, marker)

    if len(equal) > 0 and print_equal:
        if not terse:
            print("")
            print("Datasets equal on source and target:")
            print("")
        _print_dataset_list(equal, verbose)

    if len(differing) > 0 and print_differing:
        if not terse:
            print("")
            print("Datasets differing between source and target:")
            print("")
        _print_dataset_list(differing, verbose)

    if len(missing) > 0 and print_missing:
        if not terse:
            print("")
            print("Datasets misssing on target:")
            print("")
        _print_dataset_list(missing, verbose)


class ParseJsonAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError("nargs not allowed")
        super(ParseJsonAction, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        parsed_values = json.loads(values)
        setattr(namespace, self.dest, parsed_values)


def main():
    # in order to have both:
    # * preformatted help text and ...
    # * automatic display of defaults
    class ArgumentDefaultsAndRawDescriptionHelpFormatter(
            argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class = ArgumentDefaultsAndRawDescriptionHelpFormatter)

    parser.add_argument(
        'base_uri', default='DATASETS', nargs='?', help='Path to local dataset repository.')
    parser.add_argument(
        '-q', '--query', default={}, dest='query',
        action=ParseJsonAction,
        help='Query for selecting datasets on lookup server.')
    parser.add_argument(
        '--compare-by', default={'uuid': True, 'name': True, 'type': True, 'creator_username': True},
        action=ParseJsonAction,
        help='Compare datasets by specific metadata.')

    parser.add_argument('--equal-only', action='store_true', dest='equal_only',
        default=False, help='Print only list of datsets found equal at source and target.')
    parser.add_argument('--missing-only', action='store_true', dest='missing_only',
        default=False, help='Print only list of datsets found missing at target.')
    parser.add_argument('--differing-only', action='store_true', dest='differing_only',
        default=False, help='Print only list of datsets found differing between source and target.')

    parser.add_argument('--terse', '-t', action='store_true', dest='terse',
        default=False, help='Machine-readible output.')
    parser.add_argument('--verbose', '-v', action='count', dest='verbose',
        default=False, help='Make this tool more verbose')
    parser.add_argument('--debug','-d', action='store_true',
        help='Make this tool print debug info')

    try:
        argcomplete.autocomplete(parser)
    except:
        pass
    args = parser.parse_args()

    loglevel = logging.ERROR

    if args.verbose > 0:
        loglevel = logging.WARN
    if args.verbose > 1:
        loglevel = logging.INFO
    if args.debug or (args.verbose > 2):
        loglevel = logging.DEBUG

    print_equal = True
    print_differing = True
    print_missing = True

    if args.equal_only:
        print_equal = True
        print_differing = False
        print_missing = False
    elif args.differing_only:
        print_equal = False
        print_differing = True
        print_missing = False
    elif args.missing_only:
        print_equal = False
        print_differing = False
        print_missing = True

    source_info = _direct_list(args.base_uri)
    target_info = _lookup_list(args.query)
    compare_dataset_lists(source_info, target_info,
                          marker=args.compare_by,
                          verbose=args.verbose, terse=args.terse,
                          print_equal=print_equal,
                          print_differing=print_differing,
                          print_missing=print_missing)


if __name__ == '__main__':
    main()
