#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
RASM snow extent compared with NSIDC

 * Data source: [Northern Hemisphere EASE-Grid Weekly Snow Cover and
   Sea Ice Extent Version 4](http://nsidc.org/data/nsidc-0046.html)
 * Comparison with RASM simulations.

Read NSIDC data

Data are stored in flat binary files with arrays of 721 columns by 721 rows
of 1­-byte unsigned integers. Each data file contains one array.
"""

import argparse
import datetime
import time as tm
import numpy as np
from netCDF4 import Dataset, date2num
import cdo
lcdo = cdo.Cdo()
import os
import sys
import re

REFERENCE_STRING = '0001-1-1 0:0:0'
REFERENCE_DATE = 10101                         # i.e. REFERENCE_STRING
REFERENCE_TIME = 0                             # i.e. REFERENCE_STRING
TIMEUNITS = 'days since ' + REFERENCE_STRING   # do not change (MUST BE DAYS)!
CALENDAR = 'standard'

tmpeasenc = 'easebin2nc_tmp.nc'
ni = 721
nj = 721


def main():

    latfile, lonfile, maskfile, \
        domainfile, binfiles, outfile = process_command_line()

    # parse dates
    startdates, enddates, middates = parse_dates(binfiles)

    sortinds = argsort(startdates)

    startdates = [startdates[i] for i in sortinds]
    enddates = [enddates[i] for i in sortinds]
    middates = [middates[i] for i in sortinds]
    binfiles = [binfiles[i] for i in sortinds]

    # get ord times
    times = date2num(middates, TIMEUNITS, calendar=CALENDAR)
    time_bounds = np.empty((len(times), 2))
    time_bounds[:, 0] = date2num(startdates, TIMEUNITS, calendar=CALENDAR)
    time_bounds[:, 1] = date2num(enddates, TIMEUNITS, calendar=CALENDAR)

    # read the bin files
    snowice = read_bin_files(binfiles)

    # setup temporary netcdf file
    write_nc(latfile=latfile,
             lonfile=lonfile,
             maskfile=maskfile,
             domainfile=domainfile,
             snowice=snowice,
             times=times,
             time_bounds=time_bounds)

    if domainfile:
        remap(domainfile, tmpeasenc, outfile)

    # clean up temporary files
    os.remove(tmpeasenc)

    return


def read_bin_files(binfiles):
    """ read akk bin files and put into a single array """
    snowice = np.empty((len(binfiles, ni, nj)))
    # Loop over input files
    for i, infile in enumerate(binfiles):
        temp = np.fromfile(infile, dtype=np.uint8)
        temp.resize([ni, nj])

        temp[temp == 5] = 1
        temp[temp != 1] = temp[0, 0]

        snowice[i, :, :] = temp[:, :]
    return snowice


def parse_dates(files):
    """ Parse dates from filename,
    i.e. EASE2_N25km.snowice.19890522-19890528.v04.bi
    """
    startdates = []
    enddates = []
    middates = []

    for fname in files:

        d0, d1 = re.findall(r'\d{8}', fname)
        t0 = datetime.datetime.strptime(d0, '%Y%m%d')
        t1 = datetime.datetime.strptime(d1, '%Y%m%d')
        tave = (t1-t0)/2 + t0

        startdates.append(t0)
        enddates.append(t1)
        middates.append(tave)

    return startdates, enddates, middates


def argsort(seq):
    """ Equivalent to numpy argsort """
    #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    #by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)


def write_nc(latfile="",
             lonfile="",
             maskfile="",
             domainfile="",
             snowice="",
             times="",
             time_bounds=""):
    """ Write the netcdf file (before remap) """

    infile = lonfile
    with open(infile, 'rb') as f:
        rawdata = f.read(ni * nj * 4)

    lon = np.frombuffer(rawdata, dtype=np.int32)
    lon.resize([ni, nj])
    lon = lon/100000.

    infile = latfile
    with open(infile, 'rb') as f:
        rawdata = f.read(ni * nj * 4)

    lat = np.frombuffer(rawdata, dtype=np.int32)
    lat.resize([ni, nj])
    lat = lat/100000.

    infile = maskfile
    with open(infile, 'rb') as f:
        rawdata = f.read(ni * nj)

    mask = np.frombuffer(rawdata, dtype=np.uint8)
    mask.resize([ni, nj])

    # Map Qc'ed values to their final categories. Then map everything that is
    # not snow to missing (that way it will not try to interpolate bewteen
    # categories I hope)

    f = Dataset(tmpeasenc, 'w', format='NETCDF4')
    f.createDimension('ni', ni)
    f.createDimension('nj', nj)
    f.createDimension('time', None)
    f.createDimension('nv', 2)

    time = f.createVariable('xc', 'f8', ('time',))
    time.long_name = "time"
    time.units = TIMEUNITS
    time.bounds = "time_bnds"
    time[:] = times

    time_bnds = f.createVariable('time_bnds', 'f8', ('time', 'nv',))
    time_bnds[:, :] = time_bounds

    x = f.createVariable('xc', 'f8', ('nj', 'ni',), fill_value=lon[0, 0])
    x.long_name = 'longitude of grid cell center'
    x.units = 'degrees_east'
    x[:, :] = lon

    y = f.createVariable('yc', 'f8', ('nj', 'ni',), fill_value=lat[0, 0])
    y.long_name = 'latitude of grid cell center'
    y.units = 'degrees_north'
    y[:, :] = lat

    dmask = f.createVariable('mask', 'u1', ('nj', 'ni',),
                             fill_value=mask[0, 0])
    dmask.long_name = 'domain mask'
    dmask.units = '-'
    dmask.coordinates = "xc yc"
    dmask[:, :] = mask

    snow = f.createVariable('snowice', 'u1', ('time', 'nj', 'ni', ),
                            fill_value=snowice[0, 0, 0])
    snow.long_name = 'Snow ice status'
    snow.units = '-'
    snow.coordinates = "xc yc"
    snow[:, :, :] = snowice

    # write attributes for netcdf
    f.description = 'Northern Hemisphere EASE-Grid 2.0 Weekly Snow Cover and Sea Ice Extent Version 4'
    f.reference = 'Brodzik, M. and R. Armstrong. 2013. Northern Hemisphere EASE-Grid 2.0 Weekly Snow Cover and Sea Ice Extent. Version 4. Boulder, Colorado USA: NASA DAAC at the National Snow and Ice Data Center. '
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0]  # prints the name of script used

    f.close()

    return


def process_command_line():
    """ Get command line args"""

    formatter_class = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description='process ease binary files '
                                                 'into a single netcdf file',
                                     formatter_class=formatter_class)
    parser.add_argument('--latfile', required=True, type=str,
                        help='latitude binary file')
    parser.add_argument('--lonfile', required=True, type=str,
                        help='longitude binary file')
    parser.add_argument('--maskfile', required=True, type=str,
                        help='ease mask file')
    parser.add_argument('--domainfile', type=str, default=False,
                        help='domain file for cdo remapping')
    parser.add_argument('--binfiles', required=True, nargs='+',
                        help='binary input files (should have two 8digit '
                             'datestrings in them for dates)')
    parser.add_argument('--outfile', type=str, default="out.nc",
                        help='name of output file')

    args = parser.parse_args()

    return (args.latfile, args.lonfile, args.maskfile, args.domainfile,
            args.binfiles, args.outfile)


def remap(grid, ifile, outfile):
    """ """
    lcdo.remapcon(grid, input=ifile, output=outfile)


if __name__ == "__main__":
    main()
