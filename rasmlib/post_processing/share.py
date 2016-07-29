"""
share.py
"""

from __future__ import print_function
import os
import socket
from warnings import warn
from ..utils import custom_strptime

NCOFORMAT = "%Y-%m-%d %H:%M:%S"

__GARNET_OPTS__ = {}

__SPIRIT_OPTS__ = {}

__HYDRA_OPTS__ = {'no_tmp_fl': True,
                  'debug': True}

host = os.getenv('HOST', socket.gethostname())

if 'spirit' in host.lower():
    MACH_OPTS = __SPIRIT_OPTS__
elif 'garnet' in host.lower():
    MACH_OPTS = __GARNET_OPTS__
elif 'hydra' in host.lower():
    MACH_OPTS = __HYDRA_OPTS__
else:
    warn('unknown host, using standard NCO options')
    MACH_OPTS = {}


class Histfile(object):
    """History File Object

    Histfile provides a well defined structure for storing file name and date
    information.

    Attributes
    ----------
    filename : str
        Filename with full path of history file.
    fname_format : str
        File name format string.  e.g. example.cpl.ha.%%Y-%%m-%%d-%%s.nc
    filedate : datetime.datetime
        Datetime object containing the date that the filename corrsponds to.
    year : int
        year of datetime.
    month : int
        month of datetime.
    day : int
        day of datetime.
    hour : int
        hour of datetime.
    minute : int
        minute of datetime.
    second : int
        second of datetime.
    """
    def __init__(self, filename='', fname_format=''):
        self.filename = filename
        self.fname_format = fname_format

        path, filename = os.path.split(filename)
        self.filedate = custom_strptime(filename, fname_format)

        self.year = self.filedate.year
        self.month = self.filedate.month
        self.day = self.filedate.day
        self.hour = self.filedate.hour
        self.minute = self.filedate.minute
        self.second = self.filedate.second

    def __str__(self):
        return "History File: {0}".format(self.filename)

    def __repr__(self):
        return "History File: {0}".format(self.filename)
