"""
module to support calculating climatolical statistics
"""

import xarray as xr
import numpy as np
from ..calendar import get_dpm


def season_mean(ds, calendar='standard'):
    """
    Calculate the seasonal mean from an xarray Dataset of monthly means. Weight
    the means by the number of days in each month.

    Parameters
    ----------
    ds : xarray.Dataset
        Monthly frequency Pandas DatetimeIndex.
    calendar : {'standard', 'gregorian', 'proleptic_gregorian', 'julian'}
        netCDF calendar with leap years.

    Returns
    ----------
    seasonal_mean : xarray.Dataset
        Dataset representing the seasonal mean of `ds`.
    """
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xr.DataArray(get_dpm(ds.time.to_index(),
                                  calendar=calendar),
                                  coords=[ds.time],
                                  name='month_length')
    # Calculate the weights by grouping by 'time.season'
    weights = (month_length.groupby('time.season') /
               month_length.groupby('time.season').sum())

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby('time.season').sum().values,
                               np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby('time.season').sum(dim='time')


def annual_mean(ds, calendar='standard'):
    """
    Calculate the annual mean from an xarray Dataset of monthly means. Weight
    the means by the number of days in each month.

    Parameters
    ----------
    ds : xarray.Dataset
        Monthly frequency Pandas DatetimeIndex.
    calendar : {'standard', 'gregorian', 'proleptic_gregorian', 'julian'}
        netCDF calendar with leap years.

    Returns
    ----------
    seasonal_mean : xarray.Dataset
        Dataset representing the annual mean of `ds`.
    """
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xr.DataArray(get_dpm(ds.time.to_index(),
                                  calendar=calendar),
                                  coords=[ds.time],
                                  name='month_length')
    # Calculate the weights by grouping by 'time.season'
    weights = month_length / month_length.sum()

    # Calculate the weighted average
    return (ds * weights).sum(dim='time')
