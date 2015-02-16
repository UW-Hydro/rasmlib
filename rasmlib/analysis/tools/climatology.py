"""
xray extensions for calculating climatological means
"""
import xray

seasons = ('DJF', 'MAM', 'JJA', 'SON')

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28.25, 31, 30, 31, 30, 31, 31, 30,
                               31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

seas = {12: 'DJF', 1: 'DJF', 2: 'DJF',
        3: 'MAM', 4: 'MAM', 5: 'MAM',
        6: 'JJA', 7: 'JJA', 8: 'JJA',
        9: 'SON', 10: 'SON', 11: 'SON'}


def xray_monthly_means(*objects):
    """
    Given any number of xray.Datasets objects, returns new objects with
    averaged by month of year.
    """
    monthly_means = []
    aligned = xray.align(objects, join='inner', copy=True)

    for ds in aligned:
        monthly_means.append(ds.groupby('time.month').mean())

    return monthly_means


def xray_seasonal_means(*objects):
    """
    Given any number of xray.Datasets objects, returns new objects with
    averaged by season.
    """
    seasonal_means = []

    aligned = xray.align(objects, join='inner', copy=True)

    for ds in aligned:
        seasonal_means.append(ds.groupby('time.season').apply(weighted_mean))

    return seasonal_means


def xray_annual_means(*objects):
    """
    Given any number of xray.Datasets objects, returns new objects with
    averaged by year.
    """
    annual_means = []

    aligned = xray.align(objects, join='inner', copy=True)

    for ds in aligned:
        annual_means.append(ds.groupby('time.year').apply(weighted_mean))

    return annual_means


def xray_water_year_means(*objects):
    """
    Given any number of xray.Datasets objects, returns new objects with
    averaged by water year.
    """
    water_year_means = []

    aligned = xray.align(objects, join='inner', copy=True)

    for ds in aligned:
        water_year_means.append(ds.groupby('time.water_year').apply(weighted_mean))

    return water_year_means


def xray_full_means(*objects):
    """
    Given any number of xray.Datasets objects, returns new objects with
    averaged over the full time period.
    """
    full_means = []

    aligned = xray.align(objects, join='inner', copy=True)

    for ds in aligned:
        full_means.append(ds.mean())

    return full_means


def weighted_mean():
    """return a time weighted mean"""
    pass

def make_new_monthly_index(old_index):
    ns = 1e-9
    t0 = old_index.values[0]
    if type(t0) == tuple:
        d1 = t0[1]
    else:
        d1 = datetime.datetime.utcfromtimestamp(t0.astype(int)*ns)
    return pd.date_range(start=d1.strftime("%Y-%m-01"), periods=len(old_index), freq='M')

def monseasmean(data, calendar='noleap', start=None, end=None):
    """
    Calculate the seasonal means (DJF, MAM, JJA, SON)
    
    inputs:     data (xray Dataset or netCDF with time axis)
    returns:    seasonal means (xray Dataset with 4 time values)
    """
       
    # Open file if not already an xray Dataset
    if isinstance(data, (str, unicode)) and os.path.isfile(data):
        data = xray.open_dataset(data)
    elif isinstance(data, xray.dataset.Dataset):
        pass
    else:
        raise TypeError('data must be either a filename or a xray.Dataset, got %s' %type(data))

    # determine time slice to average over
    time = make_new_monthly_index(data['time'])
    a, b = time.slice_locs(start, end)
    time_inds = np.arange(a, b)
    time = time[a:b][:]
    ntime = b-a

    # create weights for time slice
    seas_inds = OrderedDict()
    weights = np.empty(ntime, dtype=np.float)
    ds_seas = np.chararray(ntime, itemsize=3)
    for i, m in enumerate(time.month):
        weights[i] = dpm[calendar][m]
        ds_seas[i] = seas[m]
    for s in seasons:
        inds = np.squeeze(np.nonzero(ds_seas==s))
        weights /= weights[inds].sum()
        seas_inds[s] = inds

    variables = OrderedDict()
    data = data.isel(time=time_inds)
    for name, da in data.noncoordinates.iteritems():
        if 'time' in da.coordinates:
            new = []
            if da.ndim == 3:
                da = weights[:, np.newaxis, np.newaxis] * da
            elif da.ndim == 4:
                da = weights[:, np.newaxis, np.newaxis, np.newaxis] * da
            elif da.ndim == 2:
                da = weights[:, np.newaxis] * da
            for i, s in enumerate(seasons):
                new.append(da.isel(time=seas_inds[s]).sum(dimension='time',
                           keep_attrs=True))
            variables[name] = xray.DataArray.concat(new, dimension='season')
    final = xray.Dataset(variables=variables,
                         attributes={'history': 'created by monseasmean'})
    final['season'].values = seasons
    
    return final
