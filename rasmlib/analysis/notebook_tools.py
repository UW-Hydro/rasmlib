
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import xray
import numpy as np

seasons = ('DJF', 'MAM', 'JJA', 'SON')

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30,
                               31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}


def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
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
    return a array of days per month corresponding to the months provided in
    `time`
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if leap_year(year, calendar=calendar):
            month_length[i] += 1
    return month_length


def season_mean(ds, calendar='standard'):
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xray.DataArray(get_dpm(ds.time.to_index(),
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


# Wrap it into a simple function
def annual_mean(ds, calendar='standard'):
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xray.DataArray(get_dpm(ds.time.to_index(),
                                  calendar=calendar),
                                  coords=[ds.time],
                                  name='month_length')
    # Calculate the weights by grouping by 'time.season'
    weights = month_length / month_length.sum()

    # Calculate the weighted average
    return (ds * weights).sum(dim='time')


def sub_plot_pcolor(lons, lats, data, title=None, cmap='Spectral_r',
                    vmin=None, vmax=None, cbar=True, cbar_location='bottom',
                    units=None, projection=None):

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    if projection is None:
        projection = {'urcrnrlat': 27.511827255753555,
                      'urcrnrlon': 16.90845094934209,
                      'llcrnrlat': 16.534986367884521,
                      'llcrnrlon': 189.2229322311162,
                      'projection': 'lcc',
                      'rsphere': 6371200.0,
                      'lon_0': -114,
                      'lat_0': 90}

    m = Basemap(**projection)
    xi, yi = m(np.squeeze(lons), np.squeeze(lats))
    sp = m.pcolormesh(xi, yi, np.squeeze(data),
                      vmin=vmin, vmax=vmax, cmap=cmap)
    m.drawparallels(np.arange(-80., 81., 20.))
    m.drawmeridians(np.arange(-180., 181., 20.))
    m.drawcoastlines(color='k', linewidth=0.25)

    if title:
        plt.title(title, size=13)
    if cbar:
        cbar = m.colorbar(location=cbar_location)
    cbar.set_label(units)

    return sp


def cmap_discretize(cmap, n=10):

    cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., n), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n + 1)]

    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d" % n,
                                              cdict, 1024)


def plot4(lons, lats, pannels, variable, titles=None,
          vmin=None, vmax=None, units=None, cmap=None,
          amin=None, amax=None, amap=None,
          outdir=None, y0=None, y1=None, mask=None, relative=True, run=''):
    """
    Make 4 pannel plots for absolute values and anomaly fields
    lons - longitudes (2d array)
    lats - latitudes (2d array)
    pannels - dictionary of data pannels[pannel_number][variable][season, y, x]
    titles - optional names for each pannel
    vmin -
    vmax -
    """

    if titles is None:
        titles = pannels.keys()

    for i, season in enumerate(seasons):
        if y0 is not None and y1 is not None:
            title = '%s-%s %s %s' % (y0, y1, season, variable)
        else:
            title = '%s-%s' % (season, variable)
        fig = plt.figure(figsize=(10, 9))
        for j in range(1, 5):
            fig.add_subplot(2, 2, j)
            sub_plot_pcolor(lons, lats,
                            np.ma.masked_where(mask == 0,
                                               pannels[j][variable][i].values),
                            title=titles[j], vmin=vmin, vmax=vmax,
                            units=units, cmap=cmap)

        fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
        fig.tight_layout()

        if outdir:
            plt.savefig(os.path.join(outdir,
                                     '%s-%s-%s.png' % (variable, season, run)),
                        dpi=300, format='png', bbox_inches='tight')

        fig = plt.figure(figsize=(10, 9))
        for j in range(1, 5):
            fig.add_subplot(2, 2, j)
            if j == 4:
                sub_plot_pcolor(lons, lats,
                                np.ma.masked_where(mask == 0,
                                                   pannels[j][variable][i].values),
                                title=titles[j], vmin=vmin, vmax=vmax,
                                units=units, cmap=cmap)
            else:
                diff = (pannels[j][variable][i].values -
                        pannels[4][variable][i].values)
                sub_plot_pcolor(lons, lats,
                                np.ma.masked_where(mask == 0, diff),
                                title='%s-%s' % (titles[j], titles[4]),
                                vmin=amin, vmax=amax, units=units,
                                cmap=amap)
        fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
        fig.tight_layout()

        if outdir:
            plt.savefig(os.path.join(outdir,
                                     '%s-anom-%s-%s.png' % (variable, season,
                                                            run)),
                        dpi=300, format='png', bbox_inches='tight')

        if relative:
            fig = plt.figure(figsize=(10, 9))
            for j in range(1, 5):
                fig.add_subplot(2, 2, j)
                if j == 4:
                    sub_plot_pcolor(lons, lats,
                                    np.ma.masked_where(mask == 0,
                                                       pannels[j][variable][i].values),
                                    title=titles[j], vmin=vmin, vmax=vmax,
                                    units=units, cmap=cmap)
                else:
                    diff = (100 *
                            (pannels[j][variable][i].values -
                             pannels[4][variable][i].values) /
                            pannels[4][variable][i].values)
                    sub_plot_pcolor(lons, lats,
                                    np.ma.masked_where(mask == 0, diff),
                                    title='(%s-%s)/%s' % (titles[j], titles[4],
                                                          titles[4]),
                                    vmin=-100, vmax=100, units='%',
                                    cmap=amap)
            fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
            fig.tight_layout()

            if outdir:
                plt.savefig(os.path.join(outdir,
                            '%s-relative-anom-%s-%s.png' % (variable, season,
                                                            run)),
                            dpi=300, format='png', bbox_inches='tight')
