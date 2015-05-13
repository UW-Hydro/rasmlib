"""calendar.py"""
import datetime
import numpy as np
from netCDF4 import num2date, date2num

HOURSPERDAY = 24
SECSPERHOUR = 3600
MINSPERHOUR = 60
SECSPERMINUTE = 60
MINSPERDAY = HOURSPERDAY * MINSPERHOUR
SECSPERDAY = HOURSPERDAY * SECSPERHOUR
MONTHSPERYEAR = 12

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30,
                               31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}
seasons = ('DJF', 'MAM', 'JJA', 'SON')


def next_month(dateobj, calendar):
    """Shift dateobj to the start of the next month.

    Parameters
    ----------
    dateobj : datetime.datetime
        Datetime object to shift.
    calendar: {'noleap', '365_day', 'standard', 'gregorian',
               'proleptic_gregorian',  'all_leap', '366_day', '360_day'}
        Supported netCDF calendar.

    Returns
    -------
    dateobj : datetime.datetime
        Datetime object at the begining of the next month.
    """
    month = dateobj.month
    year = dateobj.year
    month += 1
    if month > MONTHSPERYEAR:
        year += 1
        month = 1
    return datetime.datetime(year, month, 1)


def prev_month(dateobj, calendar):
    """Shift dateobj to the last hour of the previous month

    Parameters
    ----------
    dateobj : datetime.datetime
        Datetime object to shift.
    calendar: {'noleap', '365_day', 'standard', 'gregorian',
               'proleptic_gregorian',  'all_leap', '366_day', '360_day'}
        Supported netCDF calendar.

    Returns
    -------
    dateobj : datetime.datetime
        Datetime object at the end of the previous month.
    """
    month = dateobj.month
    year = dateobj.year
    month -= 1
    if month == 0:
        year -= 1
        month = 12
    day = dpm[calendar][month]
    return datetime.datetime(year, month, day, 23)


def next_day(dateobj, calendar):
    """Shift dateobj to the start of the next day.

    Parameters
    ----------
    dateobj : datetime.datetime
        Datetime object to shift.
    calendar: {'noleap', '365_day', 'standard', 'gregorian',
               'proleptic_gregorian',  'all_leap', '366_day', '360_day'}
        Supported netCDF calendar.

    Returns
    -------
    dateobj : datetime.datetime
        Datetime object at the begining of the next day.
    """
    day = dateobj.day
    month = dateobj.month
    year = dateobj.year
    day += 1
    if day > dpm[calendar][month]:
        month += 1
        day = 1
        if month > 12:
            year += 1
            month = 1
    return datetime.datetime(year, month, day)


def prev_day(dateobj, calendar):
    """Shift dateobj to the last second of the previous day.
    Parameters
    ----------
    dateobj : datetime.datetime
        Datetime object to shift.
    calendar: {'noleap', '365_day', 'standard', 'gregorian',
               'proleptic_gregorian',  'all_leap', '366_day', '360_day'}
        Supported netCDF calendar.

    Returns
    -------
    dateobj : datetime.datetime
        Datetime object at the end of the previous day.
    """
    day = dateobj.day
    month = dateobj.month
    year = dateobj.year
    day -= 1
    if day == 0:
        month -= 1
        if month == 0:
            year -= 1
            month = 12
        day = dpm[calendar][month]
    return datetime.datetime(year, month, day, 23, 59, 59)

dpm = {'noleap': [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [None, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31,
                               30, 31],
       'all_leap': [None, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [None, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [None, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}


def leap_year(year, calendar='standard'):
    """Determine if year is a leap year

    Parameters
    ----------
    year : scalar
    calendar : {'standard', 'gregorian', 'proleptic_gregorian', 'julian'}
        netCDF calendar with leap years.  If a calendar is provided without
        leap years, this function returns `leap=False`.

    Returns
    ----------
    leap : bool
        True if `year` is a leap year, False if not.
    """
    leap = False
    if ((calendar in ['standard', 'gregorian',
                      'proleptic_gregorian', 'julian']) and
            (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
                (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
              (year % 100 == 0) and (year % 400 != 0) and (year < 1583)):
            leap = False
    return leap


def get_dpm(time, calendar='standard'):
    """
    Calculate the number of days in each month of a datetime index of monthly
    frequency. Return a array of days per month corresponding to the months
    provided in `time`.

    Parameters
    ----------
    time : pandas.DatetimeIndex
        Monthly frequency Pandas DatetimeIndex.
    calendar : {'standard', 'gregorian', 'proleptic_gregorian', 'julian'}
        netCDF calendar with leap years.

    Returns
    ----------
    month_length : numpy.ndarray, int
        1d array of month lengths.
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if month == 2 and leap_year(year, calendar=calendar):
            month_length[i] += 1
    return month_length


def day_of_year(index, calendar, units='days since 0001-01-01'):
    nc_dates = num2date(date2num(index.to_pydatetime(), units,
                        calendar=calendar), units, calendar=calendar)
    doy = np.array([d.timetuple()[-2] for d in nc_dates])
    return doy


def to_datetime(dates):
    return [datetime(*d.timetuple()[:-3]) for d in dates]
