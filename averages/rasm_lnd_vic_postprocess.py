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
import calendar
import datetime
import glob
import os
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description='Postprocess the daily VIC '
                                     'files from RASM simulations')
    parser.add_argument('--startdate', required=True, metavar='<start date in YYYY-MM-[DD]-[SSSSS]>',
                        help='Date stamp of first VIC file from simulation')
    parser.add_argument('--enddate', required=True, metavar='<end date in YYYY-MM-[DD]-[SSSSS]>',
                        help='Date stamp of last VIC file from simulation')
    parser.add_argument('--timestep', required=True, metavar='<hourly|daily|monthly>',
                        help='Timestep of VIC output')
    parser.add_argument('--srcdir', required=True, metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', required=True, metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--casename', required=True, metavar='<casename>',
                        help='RASM case name')
    parser.add_argument('--cal', required=True, metavar='<leap|noleap>',
                        help='calendar includes leap days or not')
    parser.add_argument('--delete', action='store_true',
                        help='delete files that have been updated or merged')
    parser.add_argument('--machine', required=True, metavar='<mac|hydro|spirit>',
                        help='computer executing script...sets nco paths')

    args = parser.parse_args()

    startdate = parse_date(args.startdate)
    enddate = parse_date(args.enddate)
    leapcal = args.cal.lower()
    timestep = args.timestep

    if leapcal not in ['leap', 'noleap']:
        sys.exit('Not a valid calendar: {}'.format(args.calendar))
    if timestep not in ['hourly', 'daily', 'monthly']:
        sys.exit('Not a valid timestep: {}'.format(args.timestep))
    if args.srcdir == args.destdir:
        sys.exit('Source directory and target directory must be different')
    if not os.path.exists(args.srcdir):
        sys.exit('Source path does not exist: {}'.format(args.srcdir))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    if leapcal == 'noleap': 
        dpm=[0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # Allows 1-based indexing
    else:
        dpm=[0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] # Allows 1-based indexing

    if args.machine=='mac':
        ncra = '/opt/local/bin/ncra'
        ncrcat = '/opt/local/bin/ncrcat'  
        ncap2 = '/opt/local/bin/ncap2'
        ncatted = '/opt/local/bin/ncatted'
        ncwa = '/opt/local/bin/ncwa'
    elif args.machine=='spirit':
        ncra = '/hafs_x86_64/ncra'
        ncrcat = '/hafs_x86_64/ncrcat'  
        ncap2 = '/hafs_x86_64/ncap2'
        ncatted = '/hafs_x86_64/ncatted'
        ncwa = '/hafs_x86_64/ncwa'
    else:
        ncra = '/usr/bin/ncra'
        ncrcat = '/usr/bin/ncrcat'  
        ncap2 = '/usr/bin/ncap2'
        ncatted = '/usr/bin/ncatted'
        ncwa = '/usr/bin/ncwa'

    tobedeleted = [];
    # Change the time stamp of the input files and fix time axis
    filedate = startdate
    
    if timestep == 'hourly':
        td = 1./24.
    elif timestep == 'daily':
        td = 1.
    elif timestep =='monthly':
        td = (filedate-datetime.datetime(filedate.year,filedate.month-1,1)).days

    firstdate = startdate + datetime.timedelta(days = -td)
    lastdate = enddate + datetime.timedelta(days = -td)

    newfiledate = filedate + datetime.timedelta(days = -td)
    while filedate <= enddate:
        if leapcal == 'noleap' and filedate.month == 2 and filedate.day == 29:
            filedate += datetime.timedelta(days=1)
        
        if timestep == 'hourly':
            filename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}-{4:02d}.nc'.
                        format(args.casename, filedate.year, filedate.month,
                               filedate.day),filedate.hour)
            newfilename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}-{4:02d}.nc'.
                           format(args.casename, newfiledate.year,newfiledate.month, 
                                  newfiledate.day,newfiledate.hour))

        elif timestep == 'daily':
            filename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}.nc'.
                        format(args.casename, filedate.year, filedate.month,
                               filedate.day))
            newfilename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}.nc'.
                           format(args.casename, newfiledate.year,newfiledate.month, 
                                  newfiledate.day))

        elif timestep == 'monthly':
            filename = ('{0}.vic.ha.{1:04d}-{2:02d}.nc'.
                        format(args.casename, filedate.year, filedate.month))
            newfilename = ('{0}.vic.ha.{1:04d}-{2:02d}.mean.nc'.
                           format(args.casename, newfiledate.year,newfiledate.month))
            daysinmonth = dpm[newfiledate.month]

            if leapcal == 'noleap' and filedate.month == 3 and td == 29:
                td = (daysinmonth-1/2.) +1.
            else:
                td = daysinmonth/2.+1. # this puts the output timestep at 1/2 way through the month

    
        fargs = [ncap2, '-O', '-s', 'time=time-%f' %(365.+td)]
        print '{:40s} -> {:40s}'.format(filename, newfilename)
        fargs.append(os.path.join(args.srcdir, filename))
        fargs.append(os.path.join(args.destdir, newfilename))
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        fargs = [ncatted]
        fargs.append(['-a', 'units,time,o,c,days since 1-1-1'])
        fargs.append(['-a', 'long_name,time,o,c,time'])
        fargs.append(['-a', 'dimensions,time,o,c,1'])
        fargs.append(['-a', 'calendar,time,o,c,noleap'])
        fargs.append(['-a', 'type_preferred,time,o,c,int'])
        fargs.append(['-a', 'title,global,o,c,{}'.format(newfilename)])
        fargs.append(os.path.join(args.destdir, newfilename))
        fargs = list(flatten(fargs))
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        tobedeleted.append(os.path.join(args.srcdir, filename))
        newfiledate = filedate
        
        if timestep == 'monthly':
            nextmonth = filedate.month+1
            if nextmonth>12:
                nextmonth=1
                nextyear=filedate.year+1
            else:
                nextyear = filedate.year
            filedate = datetime.datetime(nextyear,nextmonth,
                                         filedate.day,filedate.hour)
            td = (filedate-newfiledate).days

        else:
            filedate = filedate+datetime.timedelta(td)

    if timestep == 'hourly':
        lasthour = 23
    elif timestep == 'daily':
        lasthour = 0
        lastday = calendar.monthrange(lastdate.year, lastdate.month)[1]
    elif timestep == 'monthly':
        lasthour = 0
        lastday = 1

    # Create monthly and annual average using ncra (note that only whole months
    # and years will be averaged)
    # find first timestep in timeseries
    if timestep != 'monthly':
        hour = firstdate.hour
        day = firstdate.day
        month = firstdate.month
        year = firstdate.year
        if (hour != 0 or day != 1):
            hour = 0
            day = 1
            month = firstdate.month + 1
            year = firstdate.year
            if month > 12:
                month = 1
                year += 1
        firstdate = datetime.datetime(year, month, day, hour)
    
        # find last timestep in timeseries
        hour = lastdate.hour
        day = lastdate.day
        month = lastdate.month
        year = lastdate.year
        daysinmonth = calendar.monthrange(year, month)[1]
        if hour != lasthour:
            hour = lasthour
            day -= 1
            if day < 1:
                lastdate -= datetime.timedelta(days=1)
                hour = lasthour
                day = lastdate.day
                month = lastdate.month
                year = lastdate.year
            if leapcal == 'noleap' and month == 2:
                daysinmonth = 28
                day = 28
            if day != daysinmonth:
                month -= 1
                if month < 1:
                    year -= 1
                    month = 12
            day = calendar.monthrange(year, month)[1]
            if leapcal == 'noleap' and month == 2:
                day = 28
        lastdate = datetime.datetime(year, month, day, hour)
    
    currentdate = firstdate
    firstmonth = currentdate
    
    # monthly means
    while currentdate <= lastdate:
        outfile = ('{0}.vic.ha.{1:04d}-{2:02d}.mean.nc'.
                   format(args.casename, currentdate.year, currentdate.month))
        if timestep != 'monthly':
            fileglob = ('{0}.vic.ha.{1:04d}-{2:02d}*nc'.
                         format(args.casename, currentdate.year, currentdate.month))
            filelist = sorted(glob.glob(os.path.join(args.destdir, fileglob)))
            fargs = list(flatten([ncra, '-O', filelist,
                              os.path.join(args.destdir, outfile)]))
            print '{} [...] -> {}'.format(fargs[0], fargs[-1])
            returnval = subprocess.call(fargs)
            if returnval != 0:
                sys.exit('Error executing: {}'.format(' '.join(fargs)))

        # add dpm variable to outfile
        fargs = list(flatten([ncap2, '-O', '-s', "dpm=0.0*time+{0}".format(dpm[currentdate.month]),
                              os.path.join(args.destdir, outfile), os.path.join(args.destdir, outfile)]))
        print '{} {} {}[...] -> {}'.format(fargs[0],fargs[1],fargs[2], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        
        # add dpm attributes
        fargs = [ncatted]
        fargs.append(['-a', 'units,dpm,o,c,days'])
        fargs.append(['-a', 'long_name,dpm,o,c,Days per month'])
        fargs.append([os.path.join(args.destdir, outfile)])
        fargs = list(flatten(fargs))
        print '{} {} {}[...] -> {}'.format(fargs[0],fargs[1],fargs[2], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))

        ndays = calendar.monthrange(currentdate.year, currentdate.month)[1]
        lastmonth = currentdate
        currentdate += datetime.timedelta(days = ndays)

    # time series of monthly means
    fileglob = '{0}.vic.ha.????-??.mean.nc'.format(args.casename)
    filelist = sorted(glob.glob(os.path.join(args.destdir, fileglob)))
    outfile = ('{0}.vic.ha.{1:04d}{2:02d}-{3:04d}{4:02d}.monthly.mean.nc'.
               format(args.casename, firstmonth.year, firstmonth.month,
                      lastmonth.year, lastmonth.month))
    fargs = list(flatten([ncrcat, '-O', filelist,
                          os.path.join(args.destdir, outfile)]))
    print '{} {} {}[...] -> {}'.format(fargs[0],fargs[1],fargs[2], fargs[-1])
    returnval = subprocess.call(fargs)
    if returnval != 0:
        sys.exit('Error executing: {}'.format(' '.join(fargs)))
    tobedeleted.append(filelist)

    # averaged monthly means: only complete years
    print 'monthly means now...'
    firstyear = firstdate.year
    lastyear = lastdate.year
    if firstdate.month != 1:
        firstyear += 1
    if lastdate.month != 12:
        lastyear -= 1

    for month in xrange(1,13):
        filelist = []
        for y in xrange(firstyear, lastyear+1):
            filelist.append(os.path.join(args.destdir, '{0}.vic.ha.{1:04d}-{2:02d}.mean.nc'.format(args.casename, y, month)))
        outfile = ('{0}.vic.ha.{1:04d}-{2:04d}-{3:02d}.monthly.mean.nc'.
                   format(args.casename, firstyear, lastyear, month))
        fargs = list(flatten([ncra, '-O', filelist,
                             os.path.join(args.destdir, outfile)]))
        print '{} {} {}[...] -> {}'.format(fargs[0],fargs[1],fargs[2], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
    
    # annual means: only complete years
    firstyear = firstdate.year
    lastyear = lastdate.year
    if firstdate.month != 1:
        firstyear += 1
    if lastdate.month != 12:
        lastyear -= 1
    currentyear = firstyear
    while currentyear <= lastyear:
        
        # make a time series for each year
        tempfile = os.path.join(args.destdir, ('{0}.vic.ha.{1:04d}.ts.nc'.format(args.casename, currentyear)))
        fileglob = ('{0}.vic.ha.{1:04d}-??.mean.nc'.
                    format(args.casename, currentyear))
        filelist = sorted(glob.glob(os.path.join(args.destdir, fileglob)))
        fargs = list(flatten([ncrcat, '-O', filelist, tempfile]))
        print '{} {} {}[...] -> {}'.format(fargs[0],fargs[1],fargs[2], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))

        # use ncwa to find annual mean
        outfile = os.path.join(args.destdir, ('{0}.vic.ha.{1:04d}.mean.nc'.
                               format(args.casename, currentyear)))
        fargs = list(flatten([ncwa, '-w', 'dpm', '-a', 'time', tempfile, outfile]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        
        if timestep != 'monthly':
            # annual time series: incomplete years as well
            fileglob = ('{0}.vic.ha.{1:04d}-??-??.nc'.
                        format(args.casename, currentyear))
            filelist = sorted(glob.glob(os.path.join(args.destdir, fileglob)))
            outfile = ('{0}.vic.ha.{1:04d}.{2}.nc'.
                       format(args.casename, currentyear, timestep))
            fargs = list(flatten([ncrcat, '-O', filelist,
                                  os.path.join(args.destdir, outfile)]))
            print '{} [...] -> {}'.format(fargs[0], fargs[-1])
            returnval = subprocess.call(fargs)
            if returnval != 0:
                sys.exit('Error executing: {}'.format(' '.join(fargs)))
            tobedeleted.append(filelist)
        
        currentyear += 1

    # make seasonal averages from monthly means
    firstyear = firstdate.year
    lastyear = lastdate.year
    if firstdate.month != 1:
        firstyear += 1
    if lastdate.month != 12:
        lastyear -= 1    

    seasons = {
    'DJF':[12, 1, 2],
    'MAM':[3, 4, 5],
    'JJA' :[6, 7, 8],
    'SON':[9, 10, 11]
    }

    for seas, mths in seasons.iteritems():
        filelist = []
        for month in mths:
            filelist.append(os.path.join(args.destdir,'{0}.vic.ha.{1:04d}-{2:04d}-{3:02d}.monthly.mean.nc'.
                   format(args.casename, firstyear, lastyear, month)))
        print filelist
        # put all months in a temp file
        tempfile = os.path.join(args.destdir,'{}.nc'.format(seas))
        fargs = list(flatten([ncrcat, '-O', filelist, tempfile]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        # get a weighted average
        outfile = os.path.join(args.destdir, ('{0}.vic.ha.{1:04d}-{2:04d}.{3}.mean.nc'.
                   format(args.casename, firstyear, lastyear, seas)))
        fargs = list(flatten([ncwa, '-w', 'dpm', '-a', 'time', tempfile, outfile]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        tobedeleted.append(tempfile)

    # cleanup
    if args.delete == True:
        print 'Cleanup: Deleting original and merged files'
        tobedeleted = list(flatten(tobedeleted))
        for file in tobedeleted:
            try:
                os.remove(file)
            except OSError as e:
                print 'Error removing file: {}'.format(e)
                pass

def flatten(*args):
    for x in args:
        if hasattr(x, '__iter__'):
            for y in flatten(*x):
                yield y
        else:
            yield x

def parse_date(datestr, sep='-'):
    try:
        temp = [int(x) for x in datestr.split(sep)]
        if len(temp)==3:
            pass
        elif len(temp)==2:
            temp.append(1)
        elif len(temp)==4:
            temp[-1] /= 3600
        else:
            raise ValueError('Length of Datestring is not between 2 and 4')
        date = datetime.datetime(*temp)

    except (ValueError, TypeError):
        sys.exit('Not a valid date: {}'.format(datestr))
    except:
        sys.exit('Unexpected error: {}'.format(sys.exc_info()[0]))
    return date

if __name__ == "__main__":
    main()
