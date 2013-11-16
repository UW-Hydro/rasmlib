#!/usr/bin/env python
"""
cdomeans.py

This script takes a single input file (timestep monthly or less)
and produces a set of climatological means

It uses cdo to calculate the means and ncks to split the temporary files

"""
import argparse
import os
from dateutil.parser import parse
from share import call_nco, MONTHSPERYEAR
from cdo import Cdo
cdo = Cdo()


def main():

    timeseries, prefix, destdir, start, end, monthly_inputs = process_command_line()

    # parse start and end
    startdate = parse(start)
    enddate = parse(end)

    outfile = os.path.join(destdir, "{}.{}-{}.nc".format(prefix, start.replace("-", ""), end.replace("-", "")))
    # shorten the file
    # Get average using ncra
    fargs = ['ncks', '-O', '-D', '0']
    fargs.append(['-d'])
    fargs.append(['time,%s,%s' %(start, end)])
    fargs.append([timeseries])
    fargs.append([outfile])
    call_nco(fargs)
    timeseries = outfile

    if not monthly_inputs:
        # Use CDO to get climate means
        # seasonal means
        seasmean = os.path.join(destdir,
                                "{}.{}-{}.{}.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), "seasmean"))
        cdo.seasmean(input=timeseries, output=seasmean)
        # Yearly Means
        yearmean = os.path.join(destdir,
                                "{}.{}-{}.{}.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), "yearmean"))
        cdo.yearmean(input=timeseries, output=yearmean)
        # Monthly Means
        ymonmean = os.path.join(destdir,
                                "{}.{}-{}.{}.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), "ymonmean"))
        cdo.ymonmean(input=timeseries, output=ymonmean)

        # Now split these files
        seasons = ['DJF', 'MAM', 'JJA', 'SON']
        print('splitting seasonal means now')
        for i, seas in enumerate(seasons):
            outfile = "{}.{}-{}.{}.mean.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), seas)
            fargs = ['ncks', '-O', '-d', 'time,%i' %(i)]
            fargs.append([seasmean, os.path.join(destdir, outfile)])
            call_nco(fargs)

        print('splitting yearly means now')
        for i, year in enumerate(xrange(startdate.year, enddate.year+1)):
            outfile = "{}.{}-{}.{}.mean.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), year)
            fargs = ['ncks', '-O', '-d', 'time,%i' %(i)]
            fargs.append([yearmean, os.path.join(destdir, outfile)])
            call_nco(fargs)

        print('splitting monthly means now')
        for i in xrange(MONTHSPERYEAR):
            month = i+1
            outfile = "{0}.{1}-{2}.{3:02d}.mean.nc".format(prefix, start.replace("-", ""), end.replace("-", ""), month)
            fargs = ['ncks', '-O', '-d', 'time,%i' %(i)]
            fargs.append([ymonmean, os.path.join(destdir, outfile)])
            call_nco(fargs)
    else:
        raise ValueError('not ready for monthly inputs')

    return

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
    parser.add_argument('timeseries', metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--prefix', help='string to prepend output files')
    parser.add_argument('--startdate', metavar='<YYYY-MM-DD>',
                        help='start date for analysis')
    parser.add_argument('--enddate', metavar='<YYYY-MM-DD>',
                        help='end date for analysis')
    parser.add_argument('--monthly_inputs', action='store_true',
                        help='indicates that data is monthly data and means should be weighted accordingly.')
    args = parser.parse_args()

    if not os.path.isfile(args.timeseries):
        raise ValueError('Timeseries file does not exist: {}'.format(args.timeseries))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    return args.timeseries, args.prefix, args.destdir, args.startdate, args.enddate, args.monthly_inputs

if __name__ == "__main__":
    main()