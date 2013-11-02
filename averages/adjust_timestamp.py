#!/usr/bin/env python
"""Generic post processor for dealing with model output with wrong timestamps

The script does the following to a series of netCDF files:
 * Fixes the time stamp of the files
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps: http://nco.sourceforge.net
"""

import argparse
import datetime
import re
import os
from share import call_nco, flatten, HOURSPERDAY, MINSPERDAY, SECSPERDAY, MONTHSPERYEAR

def main(timestep=None, nsteps=-1, fname_format=None, srcdir=None, destdir=None, calendar=None, delete=False):

    if not all([timestep, fname_format, srcdir, destdir, calendar]):
        timestep, nsteps, fname_format, srcdir, destdir, calendar, delete = process_command_line()

    if timestep == 'monthly':
        prefix = fname_format.split("%")[0]
        out_format = prefix+"%%Y-%%m.nc"
    else:
        out_format = fname_format

    # Define days in each month (Allows 1-based indexing)
    if calendar in ['noleap', '365_day']:
        days_per_month=[0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    elif calendar in ['standard', 'gregorian', 'proleptic_gregorian']:
        days_per_month=[0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    elif calendar in ['all_leap', '366_day']:
        days_per_month=[0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    elif calendar=='360_day':
        days_per_month=[0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]

    # Get a list of input files
    input_files = [ f for f in os.listdir(srcdir) if os.pathisfile(os.path.join(srcdir, f)) ].sort()
    tobedeleted = []

    # Time units (list)
    time_units = get_time_units(os.path.join(srcdir, input_files[0]))

    # Determine how to deal with the time axis
    # Amount to move the time axis
    if "days" in time_units:
        unit_mult = 1
    elif "hours" in time_units:
        unit_mult = HOURSPERDAY
    elif "minutes" in time_units:
        unit_mult = MINSPERDAY
    elif "seconds" in time_units:
        unit_mult = SECSPERDAY

    # Determine if we need to adjust time units (only if units are "xxx since 0000-1-1 0:0:0")
    if time_units[2] == "0000-1-1":
        time_units = "{} since 0001-01-01 0:0:0".format(time_units[0])
        print("adjusting netcdf base time units from 0000-1-1 to 0001-01-01 0:0:0")
        unit_offset = sum(days_per_month)*-1
    else:
        time_units = " ".join(time_units)
        unit_offset = 0.0

    # Iterate over files
    for filename in input_files:

        filedate = datetime.strptime(filename, fname_format)
        year = filedate.year
        month = filedate.month
        day = filedate.day
        hour = filedate.hour

        if calendar in ['standard', 'gregorian', 'proleptic_gregorian']:
            # Use python datetime calendar
            # Find new datestamp
            if timestep=='monthly':
                month += nsteps

                if (month+nsteps)<1:
                    year -= 1
                    month += MONTHSPERYEAR
                elif (month+nsteps)>MONTHSPERYEAR:
                    year += 1
                    month -= MONTHSPERYEAR

                newfiledate = datetime.datetime(year, month, day, hour)
                tdelta = newfiledate - filedate
                td = (tdelta.days + unit_offset) * unit_mult

            elif timestep=='daily':
                tdelta = datetime.timedelta(days=nsteps)
                newfiledate = filedate + tdelta
                td = (nsteps + unit_offset) * unit_mult

            elif timestep=='hourly':
                tdelta = datetime.timedelta(hours=nsteps)
                newfiledate = filedate + tdelta
                td = (nsteps + unit_offset) * unit_mult
        else:
            # Find new datestamp
            if timestep=='monthly':
                # if moving forward, get length of month now
                if nsteps > 0:
                    td = (days_per_month[month]+unit_offset)*unit_mult
                # determine new month and year

                month += nsteps

                if (month+nsteps)<1:
                    year -= 1
                    month += MONTHSPERYEAR
                elif (month+nsteps)>MONTHSPERYEAR:
                    year += 1
                    month -= MONTHSPERYEAR

                # if moving backwards, get length of month now
                td = (days_per_month[month]+unit_offset)*unit_mult

            elif timestep=='daily':
                # determine new day, month and year
                day += nsteps
                if (day+nsteps)<1:
                    month -= 1
                    if month<1:
                        year -= 1
                        month += MONTHSPERYEAR
                    day += days_per_month[month]
                elif (day+nsteps)>days_per_month[month]:
                    day -= days_per_month[month]
                    month += 1
                    if month>MONTHSPERYEAR:
                        year += 1
                        month -= MONTHSPERYEAR
                    day -= days_per_month[month]
                newfiledate = datetime.datetime(year,month, day, hour)

                # determine offset for time axis
                td = (nsteps+unit_offset)*unit_mult

            elif timestep=='hourly':
                # determine new hour, day, month and year
                hour += nsteps
                if (hour+nsteps)<1:
                    day -= 1
                    if day<1:
                        month -=1
                        if month<1:
                            year -= 1
                            month += MONTHSPERYEAR
                        day += days_per_month[month]
                    hour += HOURSPERDAY
                elif (hour+nsteps)>HOURSPERDAY:
                    day += 1
                    if day>days_per_month[month]:
                        day -= days_per_month[month]
                        month += 1
                        if month>MONTHSPERYEAR:
                            year += 1
                            month -= MONTHSPERYEAR
                        day -= days_per_month[month]
                    hour -= HOURSPERDAY

                # determine offset for time axis
                td = (nsteps+unit_offset)*unit_mult

            newfiledate = datetime.datetime(year,month, day, hour)

        newfilename = newfiledate.strftime(out_format)

        fargs = ['ncap2', '-O', '-s', 'time=time-%f' %(td)]
        print('{:40s} -> {:40s}'.format(filename, newfilename))
        fargs.append(os.path.join(srcdir, filename))
        fargs.append(os.path.join(destdir, newfilename))

        (out, error) = call_nco(fargs)

        fargs = ['ncatted']
        fargs.append(['-a', 'units,time,o,c,{}'.format(time_units)])
        fargs.append(['-a', 'long_name,time,o,c,time'])
        fargs.append(['-a', 'dimensions,time,o,c,1'])
        fargs.append(['-a', 'calendar,time,o,c,noleap'])
        fargs.append(['-a', 'type_preferred,time,o,c,int'])
        fargs.append(['-a', 'title,global,o,c,{}'.format(newfilename)])
        fargs.append(os.path.join(destdir, newfilename))

        (out, error) = call_nco(fargs)

        tobedeleted.append(os.path.join(srcdir, filename))

    # cleanup
    if delete:
        print('Cleanup: Deleting original and merged files')
        tobedeleted = list(flatten(tobedeleted))
        for fname in tobedeleted:
            try:
                os.remove(fname)
            except OSError as e:
                print('Error removing file: {}'.format(fname))
                print(e)
                pass

def get_time_units(filename):
    """parse the output of ncdump to find the time units"""

    fargs = ['ncdump', '-h', filename]
    (out, error) = call_nco(fargs)

    time_units = re.search('time:units = "(.*)" ;', out).group(1)

    return time_units.split()

def process_command_line():
    """ Get command line args"""

    description = """
Generic post processor for dealing with model output with wrong timestamps

The script does the following to a series of netCDF files:
 * Fixes the time stamp of the files
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps: http://nco.sourceforge.net
"""

    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--timestep', required=True, metavar='<hourly|daily|monthly>',
                        help='Timestep of time dimension in input files')
    parser.add_argument('--nsteps', default=-1, type=int,
                        help='Number of timesteps to adjust date, default=-1')
    parser.add_argument('--srcdir', required=True, metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', required=True, metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--calendar', required=True, metavar='<standard|gregorian|proleptic_gregorian|noleap|365_day|360_day|julian|all_leap|366_day>',
                        help='calendar includes leap days or not')
    parser.add_argument('--delete', action='store_true',
                        help='delete files that have been updated or merged')
    parser.add_argument('--format', required=True,
                        help='Python datetime string format (e.g. casename.component.ha.%%Y-%%m-%%d.nc)')

    args = parser.parse_args()
    calendar = args.cal.lower()
    timestep = args.timestep.lower()
    nsteps = args.nsteps

    if calendar not in ['standard', 'gregorian', 'proleptic_gregorian' 'noleap', '365_day', '360_day', 'julian', 'all_leap', '366_day']:
        raise ValueError('Not a valid calendar: {}'.format(args.calendar))
    if timestep not in ['hourly', 'daily', 'monthly']:
        raise ValueError('Not a valid timestep: {}'.format(args.timestep))
    if args.srcdir==args.destdir:
        raise ValueError('Source directory and target directory must be different')
    if not os.path.exists(args.srcdir):
        raise ValueError('Source path does not exist: {}'.format(args.srcdir))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)
    if timestep=='monthly' and abs(nsteps)>MONTHSPERYEAR:
        raise ValueError("|nsteps| must be less than 12 for timestep of monthly")
    if timestep=='daily' and abs(nsteps)>28:
        raise ValueError("|nsteps| must be less than 28 for timestep of daily")
    if timestep=='hourly' and abs(nsteps)>HOURSPERDAY:
        raise ValueError("|nsteps| must be less than 24 for timestep of hourly")

    return timestep, nsteps, args.format, args.srcdir, args.destdir, calendar, args.delete

if __name__=="__main__":
    main()
