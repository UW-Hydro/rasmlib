#!/usr/bin/env python
"""
means.py
"""
import os
from itertools import groupby
from share import dpm, next_month, next_day, prev_month, prev_day
from nco import Nco


def mean_monthly_diurnal_cycle(filelist, options, variables=None):
    """
    Creates a series of netcdf files for each year, month, and hour
    Note: only accepts hourly inputs
    """
    nco = Nco(overwrite=True)

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['outdir']
    if timestep != 'hourly':
        raise ValueError('mean_monthly_diurnal_cycle only \
                         accepts hourly inputs')
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
    outfiles = []
    for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
        for month, mgroup in groupby(ygroup, lambda x: x.filedate.month):
            for hour, hgroup in groupby(mgroup, lambda x: x.filedate.hour):
                filename = "{}.{}.hmmdc.{}-{}-{}.nc".format(casename, model,
                                                            year, month, hour)
                outfile = os.path.join(outdir, filename)
                nco.ncra(input=hgroup, ouput=outfile, variable=variables)
                outfiles.append(outfile)
    # ---------------------------------------------------------------- #
    return outfiles


def daily_mean_timeseries(filelist, options, variables=None):
    """Create a timeseries of daily means"""
    nco = Nco(overwrite=True)

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['outdir']
    tempdir = options['tempdir']
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
    daily_means = []
    if timestep == 'hourly':
        for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
            for month, mgroup in groupby(ygroup, lambda x: x.filedate.month):
                for day, dgroup in groupby(mgroup, lambda x: x.filedate.day):
                    filename = "{}.{}.hdm.{}-{}-{}.nc".format(casename, model,
                                                              year, month, day)
                    outfile = os.path.join(tempdir, filename)
                    nco.ncra(input=dgroup, ouput=outfile, variable=variables)
                    daily_means.append(outfile)
    else:
        daily_means = filelist
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Combine monthly means into a single timeseries
    format = '%Y%m'
    start = startdate.strftime(format)
    end = enddate.strftime(format)
    filename = "{}.{}.hdm.{}-{}.nc".format(casename, model, start, end)
    outfile = os.path.join(outdir, filename)
    nco.ncrcat(input=daily_means, ouput=outfile, variable=variables)
    # ---------------------------------------------------------------- #

    return outfile


def monthly_mean_timeseries(filelist, options, variables=None):
    """ Create a timeseries of monthly means"""
    nco = Nco(overwrite=True)

    # ---------------------------------------------------------------- #
    # Unpack needed options
    calendar = options['calendar']
    casename = options['casename']
    model = options['model']
    timestep = options['timestep']
    outdir = options['outdir']
    tempdir = options['tempdir']
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
    monthly_means = []
    if timestep != 'monthly':
        for year, ygroup in groupby(filelist, lambda x: x.filedate.year):
            for month, mgroup in groupby(ygroup, lambda x: x.filedate.month):
                filename = "{}.{}.hmm.{}-{}.nc".format(casename, model,
                                                       year, month)
                outfile = os.path.join(tempdir, filename)
                nco.ncra(input=mgroup, ouput=outfile, variable=variables)
                monthly_means.append(outfile)
    else:
        monthly_means = filelist
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Combine monthly means into a single timeseries
    format = '%Y%m%d'
    start = startdate.strftime(format)
    end = enddate.strftime(format)
    filename = "{}.{}.hmm.{}-{}.nc".format(casename, model, start, end)
    outfile = os.path.join(outdir, filename)
    nco.ncrcat(input=monthly_means, ouput=outfile, variable=variables)
    # ---------------------------------------------------------------- #

    return outfile
