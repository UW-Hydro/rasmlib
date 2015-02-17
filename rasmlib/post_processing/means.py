#!/usr/bin/env python
"""
means.py
"""
import os
import multiprocessing
from itertools import groupby
from .share import dpm, next_month, next_day, prev_month, prev_day, chunks
from .share import MACH_OPTS
from nco import Nco
nco = Nco(**MACH_OPTS)

global results_list
results_list = []

# Sort and grouper key functions
ykeyfunc = lambda x: x.filedate.year
mkeyfunc = lambda x: x.filedate.month
dkeyfunc = lambda x: x.filedate.day
hkeyfunc = lambda x: x.filedate.hour


def monthly_mean_diurnal_cycle(filelist, options, variables=None):
    """
    Creates a series of netcdf files for each year, month, and hour
    Note: only accepts hourly inputs
    """
    print('Making monthly diurnal cycle means')

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['directories']['monthly_mean_diurnal_cycle']
    numofproc = options['numofproc']
    if numofproc < 1:
        numofproc = 1

    if timestep not in ['hourly', 'hourlyi']:
        raise ValueError('mean_monthly_diurnal_cycle only accepts hourly '
                         'inputs')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find first and last timestamp
    startdate = filelist[0].filedate
    enddate = filelist[-1].filedate
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make sure start and end dates are at the "start" or "end" of them month
    # adjust if needed
    adjust_file_list = False
    if (startdate.day != 1) and (startdate.hour != 0):
        startdate = next_month(startdate, calendar)
        adjust_file_list = True
    if (enddate.day != dpm[calendar][enddate.month]) and (enddate.hour != 23):
        enddate = prev_month(enddate, calendar)
        adjust_file_list = True
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Update filelist
    if adjust_file_list:
        new_filelist = []
        for f in filelist:
            if (f.filedate >= startdate) and (f.filedate <= enddate):
                new_filelist.append(f)
        filelist = new_filelist
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create mean monthly diurnal cycle
    global results_list
    results_list = []

    pool = multiprocessing.Pool(numofproc)
    filelist.sort(key=ykeyfunc)
    for year, ygroup in groupby(filelist, ykeyfunc):
        ylist = list(ygroup)
        ylist.sort(key=mkeyfunc)
        for month, mgroup in groupby(ylist, mkeyfunc):
            mlist = list(mgroup)
            mlist.sort(key=hkeyfunc)
            for hour, hgroup in groupby(mlist, hkeyfunc):
                hlist = list(hgroup)
                hlist.sort(key=hkeyfunc)
                pool.apply_async(diurnal_cycle,
                                 callback=store_result,
                                 args=(year, month, hour, hlist, 
                                       outdir, casename, model,
                                       variables))
    pool.close()
    pool.join()

    # make sure the results_list is sorted
    results_list.sort()

    # ---------------------------------------------------------------- #
    return results_list


def diurnal_cycle(year, month, hour, hlist, outdir, casename, model,
                  variables):
    """ """
    filename = "{0}.{1}.hhmm.{2:04}-{3:02}-{4:02}.nc".format(casename, model,
                                                             year, month, hour)
    inputs = [fname.filename for fname in hlist]
    outfile = os.path.join(outdir, filename)
    nco.ncra(input=inputs, output=outfile, variable=variables, history=True)
    return outfile


def daily_mean_timeseries(filelist, options, variables=None):
    """Create a timeseries of daily means"""

    print('Making daily mean timeseries')

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['directories']['daily_mean_timeseries']
    tempdir = options['directories']['temp']
    numofproc = options['numofproc']
    if numofproc < 1:
        numofproc = 1
    if timestep == 'monthly':
        raise ValueError('daily_mean_timeseries requires a \
                         timestp of daily or less')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find first and last timestamp
    startdate = filelist[0].filedate
    enddate = filelist[-1].filedate
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make sure start and end dates are at the "start" or "end" of them month
    # adjust if needed
    adjust_file_list = False
    if timestep == 'hourly':
        if (startdate.day != 1) and (startdate.hour != 0):
            startdate = next_day(startdate, calendar)
            adjust_file_list = True
        maxday = dpm[calendar][enddate.month]
        if (enddate.day != maxday) and (enddate.hour != 23):
            enddate = prev_day(enddate, calendar)
            adjust_file_list = True
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Update filelist
    if adjust_file_list:
        new_filelist = []
        for f in filelist:
            if (f.filedate >= startdate) and (f.filedate <= enddate):
                new_filelist.append(f)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create daily means
    print('timestep is {0}'.format(timestep))
    if timestep == 'hourly':
        global results_list
        results_list = []
        pool = multiprocessing.Pool(numofproc)
        for year, ygroup in groupby(filelist, ykeyfunc):
            ylist = list(ygroup)
            ylist.sort(key=mkeyfunc)
            for month, mgroup in groupby(ylist, mkeyfunc):
                mlist = list(mgroup)
                mlist.sort(key=dkeyfunc)
                for day, dgroup in groupby(mlist, dkeyfunc):
                    dlist = list(dgroup)
                    dlist.sort(key=dkeyfunc)
                    pool.apply_async(day_mean,
                                     callback=store_result,
                                     args=(year, month, day, dlist,
                                           tempdir, casename,
                                           model, variables))
        pool.close()
        pool.join()

        # make sure the results_list is sorted
        results_list.sort()

    else:
        results_list = [fname.filename for fname in filelist]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # check if file list is too long, if it is chunk it up
    if len(results_list) > 1000:
        print('list of daily means is very long, chunking into sub files')
        chunked_list = chunks(results_list, 365)
        results_list = []
        pool = multiprocessing.Pool(numofproc)
        for i, chunk in enumerate(chunked_list):
            outfile = os.path.join(tempdir, 'chunk.{0}.nc'.format(i))
            pool.apply_async(cat_chunk,
                             args=(chunk, outfile))
            results_list.append(outfile)
        pool.close()
        pool.join()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Combine daily means into a single timeseries
    format = '%Y%m%d'
    start = startdate.strftime(format)
    end = enddate.strftime(format)
    filename = "{0}.{1}.hdm.{2}-{3}.nc".format(casename, model, start, end)
    outfile = os.path.join(outdir, filename)
    print('Daily mean timeseries file: {0}'.format(outfile))
    nco.ncrcat(input=results_list, output=outfile, variable=variables)
    # ---------------------------------------------------------------- #
    return outfile


def cat_chunk(chunk, outfile):
    nco.ncrcat(input=chunk, output=outfile, history=True)
    return


def day_mean(year, month, day, dlist, tempdir, casename, model,
             variables):
    filename = "{0}.{1}.hdm.{2:04}-{3:02}-{4:02}.nc".format(casename, model,
                                                            year, month, day)
    inputs = [fname.filename for fname in dlist]
    outfile = os.path.join(tempdir, filename)
    nco.ncra(input=inputs, output=outfile, variable=variables, history=True)
    return outfile


def monthly_mean_timeseries(filelist, options, variables=None):
    """ Create a timeseries of monthly means"""

    print('Making monthly mean timeseries')

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['directories']['monthly_mean_timeseries']
    tempdir = options['directories']['temp']
    numofproc = options['numofproc']
    if numofproc < 1:
        numofproc = 1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find first and last timestamp
    startdate = filelist[0].filedate
    enddate = filelist[-1].filedate
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make sure start and end dates are at the "start" or "end" of them month
    # adjust if needed
    adjust_file_list = False
    if timestep == 'hourly':
        if (startdate.day != 1) and (startdate.hour != 0):
            startdate = next_month(startdate, calendar)
            adjust_file_list = True
        maxday = dpm[calendar][enddate.month]
        if (enddate.day != maxday) and (enddate.hour != 23):
            enddate = prev_month(enddate, calendar)
            adjust_file_list = True
    elif timestep == 'daily':
        if (startdate.day != 1):
            startdate = next_month(startdate, calendar)
            adjust_file_list = True
        if (enddate.day != dpm[calendar][enddate.month]):
            enddate = prev_month(enddate, calendar)
            adjust_file_list = True
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Update filelist
    if adjust_file_list:
        new_filelist = []
        for f in filelist:
            if (f.filedate >= startdate) and (f.filedate <= enddate):
                new_filelist.append(f)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create monthly means
    global results_list
    results_list = []
    if timestep != 'monthly':
        pool = multiprocessing.Pool(numofproc)
        for year, ygroup in groupby(filelist, ykeyfunc):
            ylist = list(ygroup)
            ylist.sort(key=mkeyfunc)
            for month, mgroup in groupby(ylist, mkeyfunc):
                mlist = list(mgroup)
                mlist.sort(key=hkeyfunc)
                pool.apply_async(month_mean,
                                 callback=store_result,
                                 args=(year, month, mlist, tempdir,
                                       outdir, casename, model,
                                       variables))
        pool.close()
        pool.join()

        # make sure the results_list is sorted
        results_list.sort()
    else:
        results_list = [fname.filename for fname in filelist]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Combine monthly means into a single timeseries
    format = '%Y%m'
    start = startdate.strftime(format)
    end = enddate.strftime(format)
    filename = "{0}.{1}.hmm.{2}-{3}.nc".format(casename, model, start, end)
    outfile = os.path.join(outdir, filename)
    print('Monthly mean timeseries file: {0}'.format(outfile))
    nco.ncrcat(input=results_list, output=outfile, variable=variables)
    # ---------------------------------------------------------------- #
    return outfile


def month_mean(year, month, mlist, tempdir, outdir, casename, model,
               variables):
    filename = "{0}.{1}.hmm.{2:04}-{3:02}.nc".format(casename, model,
                                                     year, month)
    inputs = [fname.filename for fname in mlist]
    outfile = os.path.join(tempdir, filename)
    nco.ncra(input=inputs, output=outfile, variable=variables, history=True)
    return outfile


def store_result(result):
    global results_list
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    results_list.append(result)
# -------------------------------------------------------------------- #
