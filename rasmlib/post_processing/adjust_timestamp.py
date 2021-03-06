#!/usr/bin/env python
"""Generic post processor for dealing with model output with wrong timestamps

The script does the following to a series of netCDF files:
 * Fixes the time stamp of the files
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps:
http://nco.sourceforge.net
"""

from __future__ import print_function
import argparse
import datetime
import os
import glob
from warnings import warn
from nco import Nco
from ..calendar import dpm, HOURSPERDAY, MINSPERDAY, SECSPERDAY, MONTHSPERYEAR
from ..io import get_time_units
from ..utils import custom_strftime, clean_file
from .share import Histfile, MACH_OPTS

nco = Nco(**MACH_OPTS)


def main():
    """main wrapper for adjust timestep."""
    timestep, nsteps, fname_format,\
        srcdir, destdir, calendar, delete = _process_command_line()

    # ---------------------------------------------------------------- #
    # Get a list of input files
    globformat = fname_format.replace('%Y', '????')
    globformat.replace('%m', '??')
    globformat.replace('%d', '??')
    globformat.replace('%s', '?????')

    files = glob.glob(os.path.join(srcdir, globformat))
    filelist = [Histfile(filename=f, fname_format=fname_format) for f in files]
    filelist.sort(key=lambda r: r.filedate)

    print('Found %s files in %s', len(filelist), srcdir)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    newfilelist = adjust_timestamp(filelist,
                                   timestep=timestep,
                                   nsteps=nsteps,
                                   fname_format=fname_format,
                                   destdir=destdir,
                                   calendar=calendar)

    print('New file list: {0}'.format(newfilelist))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    if delete:
        for filename in filelist:
            clean_file(filename)
    # ---------------------------------------------------------------- #

    return


def adjust_timestamp(filelist,
                     timestep='',
                     nsteps=0,
                     fname_format='',
                     destdir='./',
                     calendar='standard'):
    """Adjust the timestep for all of the input files

    Parameters
    ----------
    filelist : list like.
        Iterable of `Histfile` objects.
    timestep : {'hourly', 'daily', 'monthly'}, required
        String denoting what the history file frequency was.
    nsteps : int
        Number of `timestep`s to adjust.
    fname_format : str, required
        Format string representing the filename
        (e.g. example.cpl.ha.%%Y-%%m.nc)
    destdir : str, required
        Path to where adjusted files should be written.
    calendar : {'noleap', '365_day', 'standard', 'gregorian',
               'proleptic_gregorian', 'all_leap', '366_day', '360_day'}
        Supported netCDF calendar.

    Returns
    -------
    newfilelist : list
        List of `Histfile` objects with adjusted timestamps.
    """
    print('Adjusting timestamp for {0} files...'.format(len(filelist)))
    print('Adjusting by {0} timesteps of size {1}'.format(nsteps, timestep))
    newfilelist = []

    # ---------------------------------------------------------------- #
    # Validate input parameters
    if timestep not in ['hourly', 'daily', 'monthly']:
        raise ValueError('Must provide a valid timestep.')
    if not nsteps:
        raise ValueError('Must provide a nonzero integer for nsteps.')
    if not fname_format:
        raise ValueError('Must provide a valid filename format (fname_format)')
    if not destdir:
        raise ValueError('Must provide a valid destination directory '
                         '(destdir).')
    if destdir is './':
        warn('Adjusted files will be written to current working directory.')
    if calendar not in dpm.keys():
        raise ValueError('Unknow calendar: {0}, supported calendars are: '
                         '{1}'.format(calendar, dpm.keys()))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # set out_format
    if timestep == 'monthly':
        prefix = fname_format.split("%")[0]
        out_format = prefix + "%Y-%m.nc"
    else:
        out_format = fname_format
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    days_per_month = dpm[calendar]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time units (list)
    time_units = get_time_units(filelist[0].filename).split()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine how to deal with the time axis
    # Amount to move the time axis
    if "days" in time_units:
        unit_mult = 1.
    elif "hours" in time_units:
        unit_mult = HOURSPERDAY
    elif "minutes" in time_units:
        unit_mult = MINSPERDAY
    elif "seconds" in time_units:
        unit_mult = SECSPERDAY
    else:
        unit_mult = 1.
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine if we need to adjust time units
    # (only if units are "xxx since 0000-1-1 0:0:0")
    if "0000-1-1" == time_units[2]:
        time_units = "{0} since 0001-01-01 00:00:00".format(time_units[0])
        print("adjusting netcdf base time units \
              from 0000-1-1 to 0001-01-01 00:00:00")
        unit_offset = -1 * sum(days_per_month[1:])
    else:
        time_units = None
        unit_offset = 0.0
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Iterate over files
    for fileobj in filelist:
        filename = fileobj.filename
        filedate = fileobj.filedate
        year = fileobj.year
        month = fileobj.month
        day = fileobj.day
        hour = fileobj.hour

        if calendar in ['standard', 'gregorian', 'proleptic_gregorian']:
            # Use python datetime calendar
            # Find new datestamp
            if timestep == 'monthly':
                month += nsteps

                if month < 1:
                    year -= 1
                    month += MONTHSPERYEAR
                elif month > MONTHSPERYEAR:
                    year += 1
                    month -= MONTHSPERYEAR

                newfiledate = datetime.datetime(year, month, day, hour)
                tdelta = newfiledate - filedate
                td = (tdelta.days + unit_offset) * unit_mult

            elif timestep == 'daily':
                tdelta = datetime.timedelta(days=nsteps)
                newfiledate = filedate + tdelta
                td = (nsteps + unit_offset) * unit_mult

            elif timestep == 'hourly':
                tdelta = datetime.timedelta(hours=nsteps)
                newfiledate = filedate + tdelta
                td = (nsteps + unit_offset) * unit_mult
        else:
            # Find new datestamp
            if timestep == 'monthly':
                day = 1
                hour = 0

                # get the number of days between the old date and the new one
                days = 0
                for i in range(abs(nsteps)):
                    if nsteps > 0:
                        days += days_per_month[month]
                        month += 1
                        if month > MONTHSPERYEAR:
                            year += 1
                            month -= MONTHSPERYEAR
                    else:
                        month -= 1
                        if month < 1:
                            year -= 1
                            month += MONTHSPERYEAR
                        days -= days_per_month[month]

                td = (days + unit_offset) * unit_mult

            elif timestep == 'daily':
                # determine new day, month and year
                day += nsteps
                if day < 1:
                    month -= 1
                    if month < 1:
                        year -= 1
                        month += MONTHSPERYEAR
                    day += days_per_month[month]
                elif day > days_per_month[month]:
                    day -= days_per_month[month]
                    month += 1
                    if month > MONTHSPERYEAR:
                        year += 1
                        month -= MONTHSPERYEAR
                    day -= days_per_month[month]
                newfiledate = datetime.datetime(year, month, day, hour)

                # determine offset for time axis
                td = (nsteps + unit_offset) * unit_mult

            elif timestep == 'hourly':
                # determine new hour, day, month and year
                hour += nsteps
                if (hour) < 0:
                    day -= 1
                    if day < 1:
                        month -= 1
                        if month < 1:
                            year -= 1
                            month += MONTHSPERYEAR
                        day += days_per_month[month]
                    hour += HOURSPERDAY
                elif (hour) > HOURSPERDAY - 1:
                    day += 1
                    if day > days_per_month[month]:
                        day -= days_per_month[month]
                        month += 1
                        if month > MONTHSPERYEAR:
                            year += 1
                            month -= MONTHSPERYEAR
                        day -= days_per_month[month]
                    hour -= HOURSPERDAY

                # determine offset for time axis
                td = (nsteps + unit_offset) * unit_mult

            newfiledate = datetime.datetime(int(year), int(month),
                                            int(day), int(hour))

        newfilename = os.path.join(destdir,
                                   custom_strftime(newfiledate, out_format))

        print('{0}-->{1}'.format(filename, newfilename))
        nco.ncap2(input=filename,
                  output=newfilename,
                  script='"time=time{0:+}"'.format(td))
        options = []
        if time_units:
            options.extend(['-a', 'units,time,o,c,"{0}"'.format(time_units)])
        options.extend(['-a', 'long_name,time,o,c,time'])
        options.extend(['-a', 'dimensions,time,o,c,1'])
        options.extend(['-a', 'calendar,time,o,c,noleap'])
        options.extend(['-a', 'type_preferred,time,o,c,int'])
        options.extend(['-a', 'title,global,o,c,{0}'.format(newfilename)])

        nco.ncatted(input=newfilename, output=newfilename, options=options)

        # Pack up fileobj
        fileobj.filename = newfilename
        fileobj.filedate = newfiledate
        fileobj.year = newfiledate.year
        fileobj.month = newfiledate.month
        fileobj.day = newfiledate.day
        fileobj.hour = newfiledate.hour

        newfilelist.append(fileobj)
    # ---------------------------------------------------------------- #
    return newfilelist


def _process_command_line():
    """Get command line arguments"""

    description = """
Generic post processor for dealing with model output with wrong timestamps

The script does the following to a series of netCDF files:
 * Fixes the time stamp of the files
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps:
http://nco.sourceforge.net
"""
    formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter)
    parser.add_argument('--timestep', required=True,
                        metavar=['hourly', 'daily', 'monthly'],
                        help='Timestep of time dimension in input files')
    parser.add_argument('--nsteps', default=-1, type=int,
                        help='Number of timesteps to adjust date, default=-1')
    parser.add_argument('--srcdir', required=True,
                        metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', required=True,
                        metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--calendar', required=True,
                        choices=dpm.keys(),
                        help='calendar includes leap days or not')
    parser.add_argument('--delete', action='store_true',
                        help='delete files that have been updated or merged')
    parser.add_argument('--format', required=True,
                        help='Python datetime string format \
                        (e.g. casename.component.ha.%%Y-%%m-%%d.nc), \
                        note that "%%s" denotes 5 digit day-seconds')

    args = parser.parse_args()
    calendar = args.calendar.lower()
    timestep = args.timestep.lower()
    nsteps = args.nsteps

    if args.srcdir == args.destdir:
        raise ValueError('Source directory and target directory \
                         must be different')
    if not os.path.exists(args.srcdir):
        raise ValueError('Source path does not exist: {0}'.format(args.srcdir))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    if (timestep == 'monthly') and (abs(nsteps) > MONTHSPERYEAR):
        raise ValueError("|nsteps| must be less than 12 for \
                         timestep of monthly")

    if (timestep == 'daily') and (abs(nsteps) > 28):
        raise ValueError("|nsteps| must be less than 28 for \
                         timestep of daily")
    if (timestep == 'hourly') and (abs(nsteps) > HOURSPERDAY):
        raise ValueError("|nsteps| must be less than 24 for \
                         timestep of hourly")

    return timestep, nsteps, args.format, \
        args.srcdir, args.destdir, calendar, args.delete
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
