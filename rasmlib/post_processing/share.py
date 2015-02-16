"""
share.py
"""

from __future__ import print_function
try:
    from ConfigParser import SafeConfigParser
except ImportError:
    from configparser import SafeConfigParser
import os
import datetime
import re
import tarfile
import socket
from contextlib import closing

HOURSPERDAY = 24.
SECSPERHOUR = 3600.
MINSPERHOUR = 60.
SECSPERMINUTE = 60.
MINSPERDAY = HOURSPERDAY * MINSPERHOUR
SECSPERDAY = HOURSPERDAY * SECSPERHOUR
MONTHSPERYEAR = 12

NCOFORMAT = "%Y-%m-%d %H:%M:%S"

__GARNET_OPTS__ = {}

__SPIRIT_OPTS__ = {'debug': True}

__HYDRA_OPTS__ = {'no_tmp_fl': True,
                  'ram_all': True,
                  'omp_num_threads': 2}

host = os.getenv('HOST', socket.gethostname())

if 'spirit' in host.lower():
    MACH_OPTS = __SPIRIT_OPTS__
elif 'garnet' in host.lower():
    MACH_OPTS = __GARNET_OPTS__
elif 'hydra' in host.lower():
    MACH_OPTS = __HYDRA_OPTS__
else:
    print('unknown host, using standard NCO options')
    MACH_OPTS = {}

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30,
                               31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}


class Histfile(object):
    """History File object"""
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


def store_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    results_list.append(result)
# -------------------------------------------------------------------- #


def argsort(seq):
    """ Equivalent to numpy argsort """
    # http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    # by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)


def next_month(dateobj, calendar):
    """ Move to the start of the next month """
    month = dateobj.month
    year = dateobj.year
    month += 1
    if month > MONTHSPERYEAR:
        year += 1
        month = 1
    return datetime.datetime(year, month, 1)


def prev_month(dateobj, calendar):
    """ Move to the end of the previous month """
    month = dateobj.month
    year = dateobj.year
    month -= 1
    if month == 0:
        year -= 1
        month = 12
    day = dpm[calendar][month]
    return datetime.datetime(year, month, day, 23)


def next_day(dateobj, calendar):
    """ Move to the start of the next day """
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
    """ Move to the end of the previous day """
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


def clean_dir(directory):
    """ Clean all files in a directory but leave the directory"""
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        print("cleaning -----> {0}".format(file_path))

        clean_file(file_path)
    return


def clean_file(file_name):
    """ Delete the file"""
    try:
        if os.path.isfile(file_name):
            os.unlink(file_name)
    except Exception as e:
        print(e)
        print('Error cleaning file: {0}'.format(file_name))
    return


def make_directories(rundir, subdir_names):
    """Make directory structure from subdir_names"""
    if not os.path.exists(rundir):
        os.makedirs(rundir)
    paths = {}
    for s in subdir_names:
        paths[s] = os.path.join(rundir, s)
        if not os.path.exists(paths[s]):
            os.makedirs(paths[s])
    return paths


def read_config(config_file):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """
    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = {}
    for section in sections:
        options = config.options(section)
        dict2 = {}
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2
    return dict1


def config_type(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        else:
            try:
                if '.' in value:
                    return float(value)
                else:
                    return int(value)
            except:
                return value
    else:
        try:
            return map(float, val_list)
        except:
            return val_list


def multiple_replace(dict, text):
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)


def custom_strptime(date_string, format):
    """
    Parse a custom datestring

    %s is now a strptime Directive
    (day-seconds as a 5-digit zero padded integer)
    """
    if '%s' in format:
        # do the custom parsing
        # (for now assume that the format is case.mod.avgs.%Y-%m-%d-%s.suffix)
        d = re.split('[._-]', date_string)
        year = int(d[-5])
        month = int(d[-4])
        day = int(d[-3])

        datetimeobj = datetime.datetime(year, month, day) + \
            datetime.timedelta(seconds=int(d[-2]))

    else:
        # just use the datetime module
        datetimeobj = datetime.datetime.strptime(date_string, format)

    return datetimeobj


def custom_strftime(datetimeobj, format):
    """
    Return a custom datestring from a standard datetime object.

    %s is now a strptime Directive
    (day-seconds as a 5-digit zero padded integer)
    """
    if '%s' in format:
        datetimedict = {'%Y': "{0:04d}".format(int(datetimeobj.year)),
                        '%m': "{0:02d}".format(int(datetimeobj.month)),
                        '%d': "{0:02d}".format(int(datetimeobj.day)),
                        '%s': "{0:05d}".format(int(datetimeobj.hour *
                                               SECSPERHOUR +
                                               datetimeobj.minute *
                                               SECSPERMINUTE +
                                               datetimeobj.second))}

        date_string = multiple_replace(datetimedict, format)
    else:
        date_string = datetimeobj.strftime(format)

    return date_string


def make_tarfile(output_filename, source_dir):
    """Simple wrapper to create a compressed tar file at the end of the run"""
    with closing(tarfile.open(output_filename, "w:gz")) as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def partition(lst, n):
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))]
            for i in range(n)]


def chunks(l, n):
    """ Yield successive n-sized chunks from l."""
    return [l[i:i + n] for i in range(0, len(l), n)]
