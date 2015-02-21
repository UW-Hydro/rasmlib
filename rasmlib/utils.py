"""utils.py"""
import datetime
import re
import os
from .calendar import SECSPERHOUR, SECSPERMINUTE


def argsort(seq):
    """Returns the indices that would sort an iterable

    Equivalent to numpy argsort but works on generic Python interables

    Parameters
    ----------
    seq : list like
        Iterable to sort.

    Returns
    -------
    index_list : list
        Array of indices that sort `seq` along the specified axis.
        In other words, ``seq[index_array]`` yields a sorted `seq`.
    """
    # http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    # by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)


def clean_dir(directory):
    """Clean all files in a directory but leave the directory

    Parameters
    ----------
    directory : str
        Directory path to remove files from.
    """
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        print("cleaning -----> {0}".format(file_path))

        clean_file(file_path)
    return


def clean_file(file_name):
    """Delete the file

    Parameters
    ----------
    file_name : str
        Filename to remove.
    """
    try:
        if os.path.isfile(file_name):
            os.unlink(file_name)
    except Exception as e:
        print('Error cleaning file: {0}'.format(file_name))
        print(e)
    return


def make_directories(topdir, subdir_names):
    """Make directory structure from subdir_names

    Parameters
    ----------
    topdir : str
        Path to top directory
    subdir_names : list like
        List like object containing the names of subdirectories to create
        inside `topdir`

    Returns
    -------
    paths : dict
        Dictionary of full paths for each `subdir`.

    """
    if not os.path.exists(topdir):
        os.makedirs(topdir)
    paths = {}
    for s in subdir_names:
        paths[s] = os.path.join(topdir, s)
        if not os.path.exists(paths[s]):
            os.makedirs(paths[s])
    return paths


def multiple_replace(text, replace):
    """Replace sections of `text` with values in `replace`.

    Parameters
    text : str
        String of which sections will be replaced by `replace`.
    ----------
    replace : dict
        Dictionary {keys: vals} of strings to replace in `text`.  `keys`
        should be present in `text` and will be replaced with `vals`.

    Returns
    -------
    text : str
        String with replacements made.

    """
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(list(map(re.escape, replace.keys()))))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: replace[mo.string[mo.start():mo.end()]], text)


def custom_strptime(date_string, format):
    """
    Parse a custom datestring

    %s is now a strptime directive
    (day-seconds as a 5-digit zero padded integer)

    Parameters
    ----------
    date_string : str
        String to generate datetime object from (e.g. "2015-02-18").
    format : str
        Standard strptime format string.  %s may correspond to the number of
        seconds since midnight.

    Returns
    -------
    datetimeobj : datetime.datetime
        Parsed Datetime object.

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

    Parameters
    ----------
    datetimeobj : datetime.datetime
        Datetime object to convert to string using `format`.
    format : str
        Standard strptime format string.  %s may correspond to the number of
        seconds since midnight.

    Returns
    -------
    date_string : str
        Parsed type value
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

        date_string = multiple_replace(format, datetimedict)
    else:
        date_string = datetimeobj.strftime(format)

    return date_string


def partition(lst, n):
    """
    Partition the iterable `lst` into `n` pieces.

    Parameters
    ----------
    lst : list like
        List like object to be partioned.
    n : int
        Number of partitions to split `lst` into.

    Returns
    -------
    partitions : list
        `n` partitions of `lst`.
    """
    division = len(lst) / float(n)
    return [lst[int(round(division * i)): int(round(division * (i + 1)))]
            for i in range(n)]


def chunks(lst, n):
    """Yield successive n-sized chunks from `lst`.

    Parameters
    ----------
    lst : list like
        List like object to be chunked.
    n : int
        Number of chunks to split `lst` into.

    Returns
    -------
    chunks : list
        `n` chunks of `lst`.
    """
    return [lst[i:i + n] for i in range(0, len(lst), n)]
