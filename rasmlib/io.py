"""io.py"""
import os
from collections import OrderedDict
import xarray as xr
try:
    from ConfigParser import SafeConfigParser
except ImportError:
    from configparser import SafeConfigParser
import tarfile
from contextlib import closing
from netCDF4 import Dataset


class RasmLibIOError(Exception):
    pass


# -------------------------------------------------------------------- #
# Read the Configuration File
def read_config(config_filepath):
    """
    Return a dictionary with subdictionaries of all configFile options/values

    Parameters
    ----------
    config_filepath : str
        Filepath to configparser / ini style formatted configuration file

    Returns
    -------
    config_dict : collections.OrderedDict
        OrderedDict of parsed configuration file value

    """

    if not os.path.isfile(config_filepath):
        raise RasmLibIOError('Configuration File ({0}) does not '
                             'exist'.format(config_filepath))

    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_filepath)
    sections = config.sections()
    config_dict = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        config_dict[section] = dict2
    return config_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
def config_type(value):
    """Parse the type of the configuration file option.

    First see the value is a bool, then try float, finally return a string.

    If string, the `config_type` will attemp to expand any environment
    variables befoer returning.

    Parameters
    ----------
    value : str or list of strings
        Generic python value to be parsed.

    Returns
    -------
    parsed_value : {bool, int, float, str} or list of
        {bool, int, float, or str}
        Parsed value.
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
        elif _isint(value):
            return int(value)
        elif _isfloat(value):
            return float(value)
        else:
            return os.path.expandvars(value)
    else:
        try:
            return list(map(float, val_list))
        except:
            pass
        try:
            return list(map(int, val_list))
        except:
            return val_list
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _isfloat(x):
    """Test of value is a `float`

    Parameters
    ----------
    x : str
        Value to test if `float`.

    Returns
    -------
    isfloat : bool
        True if `x` can be cast to `float`, otherwise False.
    """
    try:
        float(x)
    except ValueError:
        return False
    else:
        return True
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _isint(x):
    """Test if value is an integer

    Parameters
    ----------
    x : str
        Value to test if `int`.

    Returns
    -------
    isint : bool
        True if `x` can be cast to `int`, otherwise False."""
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b
# -------------------------------------------------------------------- #


def make_tarfile(output_filename, source_dir, mode="w:gz"):
    """Simple wrapper to create a compressed tar file at the end of the run

    Parameters
    ----------
    output_filename : str
        Output tarfile name.
    source_dir : str
        Directory to tar.
    mode : str
        Mode to open and compress tarfile.  See `tarfile.open` for
        compression options.
    """
    with closing(tarfile.open(output_filename, mode)) as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
    return


def get_datasets(names, files, variables, analysis_vars, timestep):
    """
    Parse the files and variables namelists and load the files into xarray
    objects.
    """
    datasets = OrderedDict()

    if not any(names):
        return datasets

    for name in names:
        print('Getting data for {0}'.format(name))
        f = files.loc[(files['NAME'] == name) &
                      (files['TIMESTEP'] == timestep)]
        if len(f) > 1:
            raise RasmLibIOError('Union of NAME: {0} and TIMESTEP: {1} '
                                 'returned too many rows ({2}).\n'
                                 '{3}'.format(name, timestep, len(f), f))
        elif len(f) < 1:
            raise RasmLibIOError('Union of NAME: {0} and TIMESTEP: {1} '
                                 'returned no rows.\n'
                                 '{2}'.format(name, timestep, f))
        file_path = f['FILE_PATH'].values[0]
        dataset_class = f['DATASET_CLASS'].values[0]

        # read the dataset
        ds = xr.open_dataset(file_path)

        # adjust units and var names
        for var in analysis_vars:
            v = variables[variables['VARIABLE'] == var]
            units = v['UNITS_STR'].values[0]
            dsvar = v['{0}-VARNAME'.format(dataset_class)].values[0]
            mult = v['{0}-MULT'.format(dataset_class)].values[0]
            offset = v['{0}-OFFSET'.format(dataset_class)].values[0]

            # rename variable
            if dsvar != var:
                ds[dsvar] = ds[dsvar].rename(var)
            # apply multiplier
            if mult != 1.:
                ds[dsvar] *= mult
            # apply offset
            if offset != 0.:
                ds[dsvar] += offset

            # set the units attribute
            ds[dsvar].attrs['units'] = units

            # add any attributes to the dataset
            ds.attrs['analysis_name'] = name
            ds.attrs['dataset_class'] = dataset_class

        datasets[name] = ds

    return datasets


def read_domain(filepath):
    """read a CESM domain file and return an xarray dataset

    Parameters
    ----------
    filepath : str
        CESM domain filepath

    Returns
    ----------
    domain : xarray.Dataset
        Dataset with domain variables.
    """

    re = 6.37122e6

    domain = xr.open_dataset(filepath)
    domain['xc'] = domain['xc'].rename('lon')
    domain['yc'] = domain['yc'].rename('lat')
    domain['area'] *= re * re  # area in m2

    return domain


def get_time_units(filename):
    """Get the timeunits in a netCDF file: `filename`.

    Parameters
    ----------
    filename : str
        netCDF file containing `time` variable with units attribute.

    Returns
    -------
    time_units : str
        Time units attributes in `filename`.
    """
    f = Dataset(filename)
    time_units = f.variables['time'].units
    f.close()
    return time_units
