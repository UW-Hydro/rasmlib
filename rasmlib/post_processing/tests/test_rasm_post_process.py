#!/usr/local/env python
"""
test_rasmlib.py

Set to run with pytest

Usage: py.test (from rasmlib or test directory)
"""
import pytest

import os
import tempfile
import random
import netCDF4
import numpy as np
from dateutil import relativedelta
from copy import deepcopy
from means import *
from share import *
from rasm_post_process import *
from adjust_timestamp import *

stdtimeunits = 'days since 1990-01-01'
badtimeunits = "days since 0000-1-1 0:0:0"
stdcalendar = 'standard'
noleapcalendar = 'noleap'

hourlydatestrformat = 'TEST.mod.ha.%Y-%m-%d-%s.nc'
dailydatestrformat = 'TEST.mod.ha.%Y-%m-%d.nc'
monthlydatestrformat = 'TEST.mod.ha.%Y-%m.nc'

def random_datetime():
    year = random.choice(range(1950, 2000))
    month = random.choice(range(1, 12))
    day = random.choice(range(1, 28))
    hour = random.choice(range(1, 24))
    return datetime(year, month, day, hour)


@pytest.fixture
def cleandir():
    newpath = tempfile.mkdtemp()
    os.chdir(newpath)


@pytest.fixture
def tempdestdir():
    return tempfile.mkdtemp(suffix='destdir')


@pytest.fixture
def tempsrcdir():
    return tempfile.mkdtemp(suffix='srcdir')


@pytest.fixture
def random_field():
    return np.random.rand(5, 5)


@pytest.fixture
def foo_nc(random_field, tempsrcdir):
    """ a simple netCDF file with a random field"""
    filename = os.path.join(tempsrcdir, 'foo.nc')
    f = netCDF4.Dataset(filename, 'w')
    shape = random_field.shape
    dim0 = f.createDimension('dim0', shape[0])
    dim1 = f.createDimension('dim1', shape[1])
    time = f.createDimension('time', 1)
    var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
    time = f.createVariable('time', 'f8', ('time'))
    time.units = stdtimeunits
    var[:, :, :] = random_field
    time[:] = 1.0
    f.close()
    return filename


@pytest.fixture
def hourly_filelist(random_field, hourlyfilenamelist, hourlydatetimelist,
                    tempsrcdir):
    """Create a bunch of sample hourly netcdf files with real times"""
    filelist = []
    for filename, date in zip(hourlyfilenamelist, hourlydatetimelist):
        filename = os.path.join(tempsrcdir, filename)
        f = netCDF4.Dataset(filename, 'w')
        shape = random_field.shape
        dim0 = f.createDimension('dim0', shape[0])
        dim1 = f.createDimension('dim1', shape[1])
        time = f.createDimension('time', 1)
        var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
        time = f.createVariable('time', 'f8', ('time'))
        time.units = stdtimeunits
        time.calendar = noleapcalendar
        var[:, :, :] = random_field
        time[:] = netCDF4.date2num(date, stdtimeunits, calendar=noleapcalendar)
        f.close()

        filelist.append(Histfile(filename=filename,
                                 fname_format=hourlydatestrformat))
    return filelist


@pytest.fixture
def daily_filelist(random_field, dailyfilenamelist, dailydatetimelist,
                   tempsrcdir):
    """Create a bunch of sample daily netcdf files with real times"""
    filelist = []
    for filename, date in zip(dailyfilenamelist, dailydatetimelist):
        filename = os.path.join(tempsrcdir, filename)
        f = netCDF4.Dataset(filename, 'w')
        shape = random_field.shape
        dim0 = f.createDimension('dim0', shape[0])
        dim1 = f.createDimension('dim1', shape[1])
        time = f.createDimension('time', 1)
        var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
        time = f.createVariable('time', 'f8', ('time'))
        time.units = stdtimeunits
        time.calendar = noleapcalendar
        var[:, :, :] = random_field
        time[:] = netCDF4.date2num(date, stdtimeunits, calendar=noleapcalendar)
        f.close()

        filelist.append(Histfile(filename=filename,
                                 fname_format=dailydatestrformat))
    return filelist


@pytest.fixture
def monthly_filelist(random_field, monthlyfilenamelist, monthlydatetimelist,
                     tempsrcdir):
    """Create a bunch of sample monthly netcdf files with real times"""
    filelist = []
    for filename, date in zip(monthlyfilenamelist, monthlydatetimelist):
        filename = os.path.join(tempsrcdir, filename)
        f = netCDF4.Dataset(filename, 'w')
        shape = random_field.shape
        dim0 = f.createDimension('dim0', shape[0])
        dim1 = f.createDimension('dim1', shape[1])
        time = f.createDimension('time', 1)
        var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
        time = f.createVariable('time', 'f8', ('time'))
        time.units = stdtimeunits
        time.calendar = noleapcalendar
        var[:, :, :] = random_field
        time[:] = netCDF4.date2num(date, stdtimeunits, calendar=noleapcalendar)
        f.close()

        filelist.append(Histfile(filename=filename,
                                 fname_format=monthlydatestrformat))
    return filelist


@pytest.fixture
def daily_filelist_badtimeunits(random_field, dailyfilenamelist,
                                dailydatetimelist, tempsrcdir):
    """Create a bunch of sample netcdf files with real times but bad time units
    """
    filelist = []
    for filename, date in zip(dailyfilenamelist, dailydatetimelist):
        filename = os.path.join(tempsrcdir, filename)
        f = netCDF4.Dataset(filename, 'w')
        shape = random_field.shape
        dim0 = f.createDimension('dim0', shape[0])
        dim1 = f.createDimension('dim1', shape[1])
        time = f.createDimension('time', 1)
        var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
        time = f.createVariable('time', 'f8', ('time'))
        time.units = badtimeunits
        time.calendar = noleapcalendar
        var[:, :, :] = random_field
        time[:] = netCDF4.date2num(date, badtimeunits, calendar=noleapcalendar)
        f.close()

        filelist.append(Histfile(filename=filename,
                                 fname_format=dailydatestrformat))
    return filelist


@pytest.fixture
def daily_filelist_stdcal(random_field, dailyfilenamelist, dailydatetimelist,
                          tempsrcdir):
    """Create a bunch of sample netcdf files with real times and the standard
    calendar"""
    filelist = []
    for filename, date in zip(dailyfilenamelist, dailydatetimelist):
        filename = os.path.join(tempsrcdir, filename)
        f = netCDF4.Dataset(filename, 'w')
        shape = random_field.shape
        dim0 = f.createDimension('dim0', shape[0])
        dim1 = f.createDimension('dim1', shape[1])
        time = f.createDimension('time', 1)
        var = f.createVariable('random', 'f8', ('time', 'dim0', 'dim1',))
        time = f.createVariable('time', 'f8', ('time'))
        time.units = stdtimeunits
        time.calendar = stdcalendar
        var[:, :, :] = random_field
        time[:] = netCDF4.date2num(date, badtimeunits, calendar=stdcalendar)
        f.close()

        filelist.append(Histfile(filename=filename,
                                 fname_format=dailydatestrformat))
    return filelist


@pytest.fixture
def hourlydatetimelist():
    """list of hourly datetimes"""
    return [datetime.datetime(2000, 1, 1) + i*datetime.timedelta(hours=1)
            for i in range(750)]


@pytest.fixture
def dailydatetimelist():
    """list of daily datetimes"""
    return [datetime.datetime(2000, 1, 1) + i*datetime.timedelta(days=1)
            for i in range(31)]


@pytest.fixture
def monthlydatetimelist():
    """list of monthly datetimes"""
    return [datetime.datetime(2000, 1, 1)
            + relativedelta.relativedelta(months=i) for i in range(12)]


@pytest.fixture
def hourlyfilenamelist(hourlydatetimelist):
    """list of filenames using custom_strftime"""
    return [custom_strftime(d, hourlydatestrformat)
            for d in hourlydatetimelist]


@pytest.fixture
def dailyfilenamelist(dailydatetimelist):
    """list of filenames using custom_strftime"""
    return [custom_strftime(d, dailydatestrformat)
            for d in dailydatetimelist]


@pytest.fixture
def monthlyfilenamelist(monthlydatetimelist):
    """list of filenames using custom_strftime"""
    return [custom_strftime(d, monthlydatestrformat)
            for d in monthlydatetimelist]


def test_argsort():
    """check that arg sort returns expected inds"""
    mylist = [3, 4, 1, 2, 0, 4]
    inds = [4, 2, 3, 0, 1, 5]
    assert inds == argsort(mylist)
    assert [mylist[i] for i in inds] == sorted(mylist)


def test_next_month_year():
    """test that next_month aross a year boundary works"""
    d1 = datetime.datetime(2000, 12, 30)
    assert next_month(d1, noleapcalendar) == datetime.datetime(2001, 1, 1)


def test_next_month_std():
    """test stanadard application of next_month"""
    d1 = datetime.datetime(2000, 11, 5)
    assert next_month(d1, noleapcalendar) == datetime.datetime(2000, 12, 1)


def test_prev_month_year():
    """test that prev_month aross a year boundary works"""
    d1 = datetime.datetime(2000, 1, 30)
    assert prev_month(d1, noleapcalendar) == datetime.datetime(1999, 12, 31, 23)


def test_prev_month_std():
    """test stanadard application of prev_month"""
    d1 = datetime.datetime(2000, 11, 5)
    assert prev_month(d1, noleapcalendar) == datetime.datetime(2000, 10, 31, 23)


def test_next_day_year():
    """test that next_day aross a year boundary works"""
    d1 = datetime.datetime(2000, 12, 31)
    assert next_day(d1, noleapcalendar) == datetime.datetime(2001, 1, 1)


def test_next_day_std():
    """test stanadard application of next_day"""
    d1 = datetime.datetime(2000, 11, 5, 12)
    assert next_day(d1, noleapcalendar) == datetime.datetime(2000, 11, 6)


def test_prev_day_year():
    """test that prev_day aross a year boundary works"""
    d1 = datetime.datetime(2000, 1, 1, 12)
    assert prev_day(d1, noleapcalendar) == datetime.datetime(1999, 12, 31,
                                                       23, 59, 59)


def test_prev_day_std():
    """test stanadard application of prev_day"""
    d1 = datetime.datetime(2000, 1, 5, 12)
    assert prev_day(d1, noleapcalendar) == datetime.datetime(2000, 1, 4, 23, 59, 59)


def test_clean_dir():
    """test that clean dir removes tempfile"""
    testdir = tempfile.mkdtemp()
    testfile = tempfile.mkstemp(dir=testdir)[1]
    assert os.path.isfile(testfile)
    assert len(os.listdir(testdir)) == 1
    clean_dir(testdir)
    assert len(os.listdir(testdir)) == 0


def test_clean_file():
    """test that clean_file removes tempfile"""
    testfile = tempfile.mkstemp()[1]
    clean_file(testfile)
    assert os.path.isfile(testfile) == False


def test_make_directories():
    """test that make directores creates correct subdirectories"""
    testdir = tempfile.mkdtemp()
    subdirs = ['dir1', 'dir2']
    paths = make_directories(testdir, subdirs)
    for subdir, path in paths.iteritems():
        assert os.path.isdir(path)


def test_config_type_int():
    """test that config type returns int for string of 1"""
    val = config_type('1')
    assert val == 1
    assert type(val) == int


def test_config_type_float():
    """test that config type returns float for string of 1.75"""
    val = config_type('1.75')
    assert val == 1.75
    assert type(val) == float

def test_config_type_bool():
    """test that config type returns bool True, False"""
    val = config_type('True')
    assert val
    assert type(val) == bool
    val = config_type('False')
    assert val == False
    assert type(val) == bool


def test_config_type_bool():
    """test that config type returns bool None"""
    val = config_type('None')
    assert val == None


def test_custom_strptime(hourlyfilenamelist, hourlydatetimelist):
    """test standard function of custom_strptime"""
    for f, d in zip(hourlyfilenamelist, hourlydatetimelist):
        assert d == custom_strptime(f, hourlydatestrformat)


def test_custom_strftime(hourlyfilenamelist, hourlydatetimelist):
    """test standard function of custom_strftime"""
    for f, d in zip(hourlyfilenamelist, hourlydatetimelist):
        assert f == custom_strftime(d, hourlydatestrformat)


def test_histfile_init(hourlyfilenamelist):
    """test initialization of Histfile object"""
    testfiles = [Histfile(filename=f, fname_format=hourlydatestrformat)
                 for f in hourlyfilenamelist]
    assert len(testfiles) == len(hourlyfilenamelist)


def test_hisfile_dates(hourlyfilenamelist, hourlydatetimelist):
    """test initialization of Histfile object"""
    testfiles = [Histfile(filename=f, fname_format=hourlydatestrformat)
                 for f in hourlyfilenamelist]
    filedates = [h.filedate for h in testfiles]
    assert filedates == hourlydatetimelist


def test_get_time_units(foo_nc):
    """test get time units from ncdump (get_time_units function)"""
    foounits = get_time_units(foo_nc)
    assert " ".join(foounits) == stdtimeunits


def test_adjust_timestamp_h(hourly_filelist, tempsrcdir, tempdestdir):
    """adjust timestamp of hourly files"""
    oldfilelist = deepcopy(hourly_filelist)
    newfilelist = adjust_timestamp(hourly_filelist,
                                   timestep='hourly',
                                   nsteps=-1,
                                   fname_format=hourlydatestrformat,
                                   destdir=tempdestdir,
                                   calendar=noleapcalendar)
    for new, old in zip(newfilelist, oldfilelist):
        assert new.filedate == old.filedate - datetime.timedelta(hours=1)
        assert os.path.isfile(new.filename)
        assert os.path.split(old.filename)[0] == tempsrcdir
        assert os.path.split(new.filename)[0] == tempdestdir
        assert new.filename != old.filename


def test_adjust_timestamp_d(daily_filelist, tempsrcdir, tempdestdir):
    """adjust timestamp of daily files"""
    oldfilelist = deepcopy(daily_filelist)
    newfilelist = adjust_timestamp(daily_filelist,
                                   timestep='daily',
                                   nsteps=-1,
                                   fname_format=hourlydatestrformat,
                                   destdir=tempdestdir,
                                   calendar=noleapcalendar)
    for new, old in zip(newfilelist, oldfilelist):
        assert new.filedate == old.filedate - datetime.timedelta(days=1)
        assert os.path.isfile(new.filename)
        assert os.path.split(old.filename)[0] == tempsrcdir
        assert os.path.split(new.filename)[0] == tempdestdir
        assert new.filename != old.filename


def test_adjust_timestamp_m(monthly_filelist, tempsrcdir, tempdestdir):
    """adjust timestamp of monthly files"""
    oldfilelist = deepcopy(monthly_filelist)
    newfilelist = adjust_timestamp(monthly_filelist,
                                   timestep='monthly',
                                   nsteps=-1,
                                   fname_format=hourlydatestrformat,
                                   destdir=tempdestdir,
                                   calendar=noleapcalendar)
    for new, old in zip(newfilelist, oldfilelist):
        assert new.filedate == old.filedate \
            - relativedelta.relativedelta(months=1)
        assert os.path.isfile(new.filename)
        assert os.path.split(old.filename)[0] == tempsrcdir
        assert os.path.split(new.filename)[0] == tempdestdir
        assert new.filename != old.filename


def test_adjust_timestamp_stdcal(daily_filelist_stdcal, tempsrcdir,
                                 tempdestdir):
    """adjust timestamp of daily files using standard calendar"""
    oldfilelist = deepcopy(daily_filelist_stdcal)
    newfilelist = adjust_timestamp(daily_filelist_stdcal,
                                   timestep='daily',
                                   nsteps=-1,
                                   fname_format=dailydatestrformat,
                                   destdir=tempdestdir,
                                   calendar=stdcalendar)
    for new, old in zip(newfilelist, oldfilelist):
        assert new.filedate == old.filedate - datetime.timedelta(days=1)
        assert os.path.isfile(new.filename)
        assert os.path.split(old.filename)[0] == tempsrcdir
        assert os.path.split(new.filename)[0] == tempdestdir
        assert new.filename != old.filename


def test_adjust_timestamp_badtimeunits(daily_filelist_badtimeunits,
                                       tempsrcdir, tempdestdir):
    """adjust timestamp of daily files using bad time units"""
    oldfilelist = deepcopy(daily_filelist_badtimeunits)
    newfilelist = adjust_timestamp(daily_filelist_badtimeunits,
                                   timestep='daily',
                                   nsteps=-1,
                                   fname_format=dailydatestrformat,
                                   destdir=tempdestdir,
                                   calendar=stdcalendar)
    for new, old in zip(newfilelist, oldfilelist):
        new_units = get_time_units(new.filename)
        old_units = get_time_units(old.filename)
        assert new.filedate == old.filedate - datetime.timedelta(days=1)
        assert os.path.isfile(new.filename)
        assert os.path.split(old.filename)[0] == tempsrcdir
        assert os.path.split(new.filename)[0] == tempdestdir
        assert new.filename != old.filename
        assert new_units != old_units
        assert new_units[2] == "0001-01-01"

def test_read_config_file():
    """test for successful read of sample_config_file.cfg and dict return"""
    config_dict = read_config('sample_config_file.cfg')
    assert type(config_dict) == dict

