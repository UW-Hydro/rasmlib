#!/usr/bin/env python
"""
means.py
"""
import os
import multiprocessing
from itertools import groupby
from share import dpm, next_month, next_day, prev_month, prev_day, chunks
from nco import Nco
nco = Nco()

global results_list


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
    tempdir = options['directories']['temp']
    numofproc = options['numofproc']

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
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create mean monthly diurnal cycle
    results_list = []

    if numofproc > 1:
        pool = multiprocessing.Pool(numofproc)
        for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
            for month, mgroup in groupby(ygroup, lambda x: x.filedate.month):
                    pool.apply_async(diurnal_cycle,
                                     callback=store_result,
                                     args=(year, month, mgroup, tempdir,
                                           outdir, casename, model,
                                           variables))
        pool.close()
        pool.join()

        # make sure the results_list is sorted
        results_list.sort(key=lambda x: x.filedate)

    else:
        for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
            for month, mgroup in groupby(ygroup, lambda x: x.filedate.month):
                outfile = diurnal_cycle(year, month, mgroup, tempdir,
                                        outdir, casename, model,
                                        variables)
                results_list.append(outfile)
    # ---------------------------------------------------------------- #
    return results_list


def diurnal_cycle(year, month, mgroup, tempdir, outdir, casename, model,
                  variables):
    """ """
    hours = []
    for hour, hgroup in groupby(mgroup, lambda x: x.filedate.hour):
        filename = "{0}.{1}.hhmm.{2:04}-{3:02}-{4:02}.nc".format(casename,
                                                                 model, year,
                                                                 month, hour)
        inputs = [fname.filename for fname in hgroup]
        outfile = os.path.join(tempdir, filename)
        nco.ncra(input=inputs, output=outfile, variable=variables,
                 options='-h')
        hours.append(filename)

    filename = "{0}.{1}.hmmdc.{2:04}-{3:02}.nc".format(casename, model, year,
                                                       month)
    nco.ncrcat(input=hours, output=filename, variables=variables)

    return filename


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
    results_list = []
    if timestep == 'hourly':
        if numofproc > 1:
            pool = multiprocessing.Pool(numofproc)

            for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
                for month, mgroup in groupby(ygroup,
                                             lambda x: x.filedate.month):
                    for day, dgroup in groupby(mgroup,
                                               lambda x: x.filedate.day):
                        pool.apply_async(day_mean,
                                         callback=store_result,
                                         args=(year, month, mgroup, tempdir,
                                               outdir, casename, model,
                                               variables))
            pool.close()
            pool.join()

            # make sure the results_list is sorted
            results_list.sort(key=lambda x: x.filedate)

        else:
            for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
                for month, mgroup in groupby(ygroup,
                                             lambda x: x.filedate.month):
                    for day, dgroup in groupby(mgroup,
                                               lambda x: x.filedate.day):
                        outfile = day_mean(year, month, day, dgroup, tempdir,
                                           outdir, casename, model, variables)
                        results_list.append(outfile)
    else:
        results_list = [fname.filename for fname in filelist]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # check if file list is too long, if it is chunk it up
    if len(results_list) > 3650:
        print('list of daily means is very long, chunking into sub files')
    chunked_list = chunks(results_list, 365)
    results_list = []
    for i, chunk in enumerate(chunked_list):
        outfile = os.path.join(tempdir, 'chunk.{0}.nc'.format(i))
        nco.ncrcat(input=chunk, output=outfile, options='-h')
        results_list.append(outfile)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Combine monthly means into a single timeseries
    format = '%Y%m%d'
    start = startdate.strftime(format)
    end = enddate.strftime(format)
    filename = "{0}.{1}.hdm.{2}-{3}.nc".format(casename, model, start, end)
    outfile = os.path.join(outdir, filename)
    print('Daily mean timeseries file: {0}'.format(outfile))
    nco.ncrcat(input=results_list, output=outfile, variable=variables,
               omp_num_threads=numofproc)
    # ---------------------------------------------------------------- #
    return outfile


def day_mean(year, month, day, dgroup, tempdir, outdir, casename, model,
             variables):
    filename = "{0}.{1}.hdm.{2:04}-{3:02}-{4:02}.nc".format(casename, model,
                                                            year, month, day)
    inputs = [fname.filename for fname in dgroup]
    outfile = os.path.join(tempdir, filename)
    nco.ncra(input=inputs, output=outfile, variable=variables, options='-h')
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
            startdate = next_month(startdate)
            adjust_file_list = True
        maxday = dpm[calendar][enddate.month]
        if (enddate.day != maxday) and (enddate.hour != 23):
            enddate = prev_month(enddate, calendar)
            adjust_file_list = True
    elif timestep == 'daily':
        if (startdate.day != 1):
            startdate = next_month(startdate)
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
    results_list = []
    if timestep != 'monthly':
        if numofproc > 1:
            pool = multiprocessing.Pool(numofproc)
            for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
                for month, mgroup in groupby(ygroup,
                                             lambda x: x.filedate.month):
                    pool.apply_async(month_mean,
                                     callback=store_result,
                                     args=(year, month, mgroup, tempdir,
                                           outdir, casename, model,
                                           variables))
            pool.close()
            pool.join()

            # make sure the results_list is sorted
            results_list.sort(key=lambda x: x.filedate)
        else:
            for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
                for month, mgroup in groupby(ygroup,
                                             lambda x: x.filedate.month):
                    outfile = month_mean(year, month, mgroup, tempdir, outdir,
                                         casename, model, variables)
                    results_list.append(outfile)
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
    nco.ncrcat(input=results_list, output=outfile, variable=variables,
               omp_num_threads=numofproc)
    # ---------------------------------------------------------------- #
    return outfile


def month_mean(year, month, mgroup, tempdir, outdir, casename, model,
               variables):
    filename = "{0}.{1}.hmm.{2:04}-{3:02}.nc".format(casename, model,
                                                     year, month)
    inputs = [fname.filename for fname in mgroup]
    outfile = os.path.join(tempdir, filename)
    nco.ncra(input=inputs, output=outfile, variable=variables, options='-h')
    return outfile


def store_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    results_list.append(result)
# -------------------------------------------------------------------- #
