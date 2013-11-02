#!/usr/bin/env python
"""Postprocess the daily VIC files from RASM simulations

The script does the following to a series of daily RASM VIC files:
 * Fixes the time stamp of the daily VIC files
 * Calculates monthly and annual averages
 * Concatenates the daily files into a annual file that can be used for further
   processing
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps: http://nco.sourceforge.net
"""

import argparse
import re
import datetime
from dateutil.parser import parse
import os
from share import argsort, call_nco, NCOFORMAT, MONTHSPERYEAR


def main(fname_format=None, srcdir=None, destdir=None, partial_months=False):

    if not all([fname_format, srcdir, destdir]):
        fname_format, srcdir, destdir, partial_months = process_command_line()

    # Get a list of input files
    input_files = [ f for f in os.listdir(srcdir) if os.pathisfile(os.path.join(srcdir, f)) ].sort()

    # Get timesteps for all files and create a mapping from them to their file
    dates = []
    files = []
    for i, fname in enumerate(input_files):
        d = get_netcdf_datetime(os.path.join(srcdir, fname))
        dates.append(d)
        files.append([i for j in xrange(len(d))])

    # They should be sorted but you never know
    sortinds = argsort(dates)
    dates = dates[sortinds]
    files = files[sortinds]

    # Decide what to do with partial months
    if partial_months:
        start = dates[0]
        end = dates[-1]
    else:
        start = dates[0]
        end = dates[-1]
        raise

    # Set initial values
    bound0 = datetime.datetime(start.year, start.month, 1, 0)
    bound1 = next_month(bound0)

    # Determine how to format the outfiles
    prefix = fname_format.split("%")[0]
    out_format = prefix+"%%Y-%%m.nc"

    # Make the monthly averages
    while bound0 < end:

        # Get list of infiles
        dinds = dates>=bound0 and dates<bound1
        infiles = files[dinds]

        # Get outfile name
        outfile = bound0.strftime(out_format)

        # Get average using ncra
        fargs = ['ncra', '-O']
        fargs.append(['-d' 'time', bound0.strftime(NCOFORMAT), bound1.strftime(NCOFORMAT)])
        fargs.append(['--path', srcdir])
        fargs.append(infiles)
        fargs.append(outfile)
        (out, error) = call_nco(fargs)

        # Move forward
        bound0 = bound1
        bound1 = next_month(bound1)

    return

def next_month(bound):
    """ Move to the start of the next month """
    month = bound.month
    year = bound.year
    month += 1
    if month > MONTHSPERYEAR:
        year += 1
        month = 1
    return datetime.datetime(year, month, 1, 0)

def get_netcdf_datetime(filename, s=None):
    """Use ncdump -c -t to get a the time range in a netcdf file"""

    if not s:
        s = slice(None)

    fargs = ['ncdump', '-c', '-t', filename]
    (out, error) = call_nco(fargs)

    # split the output from ncdump at "time",
    # take the last list item which should be the timestamps
    # use re to find all occurances of the datestamps
    # make sure it is a list
    times = list(re.findall(r'\"(.+?)\"', out.split("time")[-1]))

    # Use the dateutil parser to return datetimes
    # Use the slice
    return [ parse(time) for time in times ][s]

def process_command_line():
    """ Get command line args"""

    description = """
Generic post processor for creating monthly mean timeseries of netCDF files

The script does the following to a series of netCDF files:
 * Creates monthly mean files from a series

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps: http://nco.sourceforge.net (requires version 4.2.1 or later)
http://nco.sourceforge.net/nco.html#Multiple-files-with-multiple-time-points
"""

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--srcdir', required=True, metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', required=True, metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--format', required=True,
                        help='Python datetime string format (e.g. casename.component.ha.%%Y-%%m-%%d.nc)')
    parser.add_argument('--partial_months', type=bool, default=False,
                        help='Allow allow partial months')

    args = parser.parse_args()

    if args.srcdir==args.destdir:
        raise ValueError('Source directory and target directory must be different')
    if not os.path.exists(args.srcdir):
        raise ValueError('Source path does not exist: {}'.format(args.srcdir))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    return args.format, args.srcdir, args.destdir, args.partial_months

if __name__ == "__main__":
    main()
