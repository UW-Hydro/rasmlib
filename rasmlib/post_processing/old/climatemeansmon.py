#!/usr/bin/env python
"""
climatemeansmon.py

This script takes a single input file (timeseries of monthly means)
and produces a set of climatological means

It uses cdo to calculate the means and ncks to split the temporary files

"""
import argparse
import os
import subprocess
from string import Template
from datetime import date
import calendar
from dateutil.parser import parse
from share import dpm
from cdo import Cdo
from nco import Nco
cdo = Cdo()
nco = Nco(debug=True, ram_all=True)

try:
    from netCDF4 import Dataset
    netCDF4_available = True
except:
    netCDF4_available = False

dimdefaults = {'min': '', 'max': '', 'stride': '', 'subcycle': ''}
timdim = Template("time,$min,$max,$stride,$subcycle")


class DateObj(object):
    def __init__(self, year, month, cal='standard'):
        self.year = year
        self.month = month
        self.calendar = cal

        # Set days in month
        if self.calendar in ['standard', 'gregorian', 'proleptic_gregorian']:
            self.days = calendar.monthrange(self.year, self.month)[1]
        elif self.calendar in dpm.keys():
            self.days = dpm[self.calendar][self.month]
        else:
            raise ValueError('Calendar {0} not in dpm '
                             'dictionary'.format(self.calendar))

        # Set Season
        if self.month in [1, 2, 12]:
            self.season = 'DJF'
        elif self.month in [3, 4, 5]:
            self.season = 'MAM'
        elif self.month in [6, 7, 8]:
            self.season = 'JJA'
        elif self.month in [9, 10, 11]:
            self.season = 'JJA'
        else:
            raise ValueError('Month {0} out of range'.format(self.month))

        # Set Water Year
        if month > 9:
            self.water_year = self.year + 1
        else:
            self.water_year = self.year

        # Set Water Year Season
        if self.month in [10, 11, 12]:
            self.season = 'OND'
        elif self.month in [1, 2, 3]:
            self.season = 'JFM'
        elif self.month in [4, 5, 6]:
            self.season = 'AMJ'
        elif self.month in [7, 8, 9]:
            self.season = 'JAS'
        else:
            raise ValueError('Month {0} out of range'.format(self.month))

    def __str__(self):
        return "DateObj({0}, {1})".format(self.year, self.month)

    def __repr__(self):
        return "DateObj({0}, {1})".format(self.year, self.month)


def main():

    timeseries, prefix,
    destdir, start, end, monthly_inputs = process_command_line()

    # parse start and end
    startdate = parse(start)
    enddate = parse(end)

    return


def add_dpm(filename, dayspermonth):
    """Add dpm variable to file"""

    if netCDF4_available:
        f = Dataset(filename, 'a')
        dpm_var = f.createVariable('dpm', 'f8', ('time', ))
        dpm_var[:] = dayspermonth
        dpm_var.long_name = "Days per month"
        dpm_var.units = "days"
        f.close()
    else:
        nctext = """netcdf dpm {
dimensions:
    time=unlimited;
variables:
    float dpm(time);
    dpm:long_name="Days per month";
    dpm:units="days";
data:
    dpm=%s;
}""" % (",".join(str(x) for x in dayspermonth))

        with open('dpm.txt', mode='w') as f:
            f.write(nctext)

        fargs = ['ncgen', '-b', '-o', 'dpm.nc', 'dpm.txt']
        returnval = subprocess.call(fargs)

        nco.ncks(input='dpm.nc',
                 output=filename,
                 append=True,
                 variable='dpm')
    return filename


def long_term_mean():
    """calculate long term mean"""
    dims = {'min': 'something',
            'max': 'something'}
    dimension = timdim.substitue(dimdefaults, dims)
    outfile = "something.nc"
    nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
             dimension=dimension)
    return outfile


def annual_means():
    """Calculate Annual Means from Monthly Data"""
    outfiles = []
    for year in years:
        dims = {'min': 'something',
                'max': 'something'}
        dimension = timdim.substitue(dimdefaults, dims)
        outfile = "something.nc"
        nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                 dimension=dimension)
        outfiles.append(outfile)
    return outfiles


def seasonal_means():
    """Calculate Seasonal Means from Monthly Data"""
    # Calculate Seasonal Means for each year
    outfiles = []
    for year in years:
        for seas in seasons:
            dims = {'min': 'something',
                    'max': 'something',
                    'stride': 12,
                    'subcycle': 3}
            dimension = timdim.substitue(dimdefaults, dims)
            outfile = "something.nc"
            nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                     dimension=dimension)
            outfiles.append(outfile)

    # Calculate Seasonal Means for all time
    for seas in seasons:
        dims = {'min': 'something like 0',
                'max': 'something like -1',
                'stride': 12,
                'subcycle': 3}
        dimension = timdim.substitue(dimdefaults, dims)
        outfile = "something.nc"
        nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                 dimension=dimension)
        outfiles.append(outfile)
    return outfiles


def water_year_means():
    """Calculate Annual Water Year Means from Monthly Data"""
    outfiles = []
    for year in years:
        dims = {'min': 'something',
                'max': 'something'}
        dimension = timdim.substitue(dimdefaults, dims)
        outfile = "something.nc"
        nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                 dimension=dimension)
        outfiles.append(outfile)
    return outfiles


def water_year_seasonal_means():
    """
    Calculate Water Year Seasonal Means from Monthly Data
    OND, JFM, AMJ, JAS
    """
    for year in years:
        outfiles = []
        for seas in wyseasons:
            dims = {'min': 'something',
                    'max': 'something',
                    'stride': 12,
                    'subcycle': 3}
            dimension = timdim.substitue(dimdefaults, dims)
            outfile = "something.nc"
            nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                     dimension=dimension)
            outfiles.append(outfile)
    return outfiles


def Nmonth_means(N):
    outfiles = []
    for year in years:
        for seas in wyseasons:
            dims = {'min': 'something',
                    'max': 'something',
                    'stride': 12,
                    'subcycle': N}
            dimension = timdim.substitue(dimdefaults, dims)
            outfile = "something.nc"
            nco.ncwa(input=filename, output=outfile, wgt_var='dpm',
                     dimension=dimension)
            outfiles.append(outfile)
    return outfiles


def process_command_line():
    """ Get command line args"""

    description = """
Generic post processor for creating monthly mean timeseries of netCDF files

The script does the following to a series of netCDF files:
 * Creates monthly mean files from a series

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps:
http://nco.sourceforge.net (requires version 4.2.1 or later)
"""

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('timeseries', metavar='<source directory>',
                        help='Input ')
    parser.add_argument('--destdir', metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--prefix', help='string to prepend output files')
    parser.add_argument('--startdate', metavar='<YYYY-MM-DD>',
                        help='start date for analysis')
    parser.add_argument('--enddate', metavar='<YYYY-MM-DD>',
                        help='end date for analysis')
    parser.add_argument('--monthly_inputs', action='store_true',
                        help='indicates that data is monthly data and means '
                             'should be weighted accordingly.')
    args = parser.parse_args()

    if not os.path.isfile(args.timeseries):
        raise ValueError('Timeseries file does not exist: '
                         '{0}'.format(args.timeseries))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    return args.timeseries, args.prefix, args.destdir, args.startdate,
    args.enddate, args.monthly_inputs

if __name__ == "__main__":
    main()
