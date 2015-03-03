"""
rasmlib plot utilities
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np
from ..calendar import seasons

projections = {}
projections['wr50a'] = {'urcrnrlat': 27.511827255753555,
                        'urcrnrlon': 16.90845094934209,
                        'llcrnrlat': 16.534986367884521,
                        'llcrnrlon': 189.2229322311162,
                        'projection': 'lcc',
                        'rsphere': 6371200.0,
                        'lon_0': -114,
                        'lat_0': 90}
default_projection = projections['wr50a']

default_map = None  # set again at end of module


class Bmap(object):
    """Wrapper container for Basemap object and map indices"""
    def __init__(self, projection=default_projection):
        self.projection = projection
        self.m = Basemap(**self.projection)
        self.inds_set = False

    def set_map_inds(self, lons, lats):
        """Set the map indices
        lons : numpy.ndarray
            Array of grid cell longitudes.
        lons : numpy.ndarray
            Array of grid cell latitudes.
        """
        self.xi, self.yi = self.m(np.squeeze(lons), np.squeeze(lats))
        self.inds_set = True


def make_bmap(projection=default_projection, lons=None, lats=None):
    """
    Make a bmap object storing the projection and basemap indices for a series
    of plots.

    Parameters
    ----------
    projection : dict
        Projection keywords to be passed to `Basemap`.
    lons : numpy.ndarray
        Array of grid cell longitudes.
    lons : numpy.ndarray
        Array of grid cell latitudes.

    Returns
    ----------
    bmap : Bmap
        Bmap object.
    """
    bmap = Bmap(projection=projection)

    if lons is not None and lats is not None:
        bmap.set_map_inds(lons, lats)
    return bmap


def sub_plot_pcolor(data,
                    title=None,
                    cmap='Spectral_r',
                    vmin=None,
                    vmax=None,
                    cbar=True,
                    cbar_location='bottom',
                    units=None,
                    cbar_extend='neither',
                    map_obj=default_map,
                    ax=None):
    """Plot data into a subplot using pcolormesh.

    Parameters
    ----------
    data : 2d numpy array
        Array of data to plot using pcolormesh.
    title : str, optional
        Title to add to subplot.
    cmap : str or colormap object, optional
        Colormap to use in pcolormesh.
    vmin : scalar, optional
        Minimum value for color range.
    vmax : scalar, optional
        Maximum value for color range.
    cbar_location : str, optional
        Location of colorbar.  Value is passed to `Basemap.colorbar()`
    cbar_extend : str, optional
        Extend colorbar ends.  Value is passed to `Basemap.colorbar()`
    units : str
        Colorbar label.
    map_obj : Bmap, optional
        `Bmap` object containing `Basemap` object as well as plot coordinates.
        Set using `make_bmap()`.
    ax : axis object, optional
        Axis object to use for subplot.  If None, `ax=plt.gca()`.
    """

    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    if ax is None:
        ax = plt.gca()

    mappable = map_obj.m.pcolormesh(map_obj.xi, map_obj.yi, np.squeeze(data),
                                    vmin=vmin, vmax=vmax, cmap=cmap, ax=ax)
    map_obj.m.drawparallels(np.arange(-80., 81., 20.))
    map_obj.m.drawmeridians(np.arange(-180., 181., 20.))
    map_obj.m.drawcoastlines(color='k', linewidth=0.25)

    if title is not None:
        plt.title(title, size=13)
    if cbar:
        cbar = map_obj.m.colorbar(location=cbar_location, extend=cbar_extend,
                                  mappable=mappable)
        cbar.set_label(units)

    return


def cmap_discretize(cmap, n_colors=10):
    """Return discretized colormap.

    Parameters
    ----------
    cmap : str or colormap object
        Colormap to discretize.
    n_colors : int
        Number of discrete colors to divide `cmap` into.

    Returns
    ----------
    disc_cmap : LinearSegmentedColormap
        Discretized colormap.
    """
    try:
        cmap = cm.get_cmap(cmap)
    except:
        cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n_colors + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n_colors + 1)]

    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d" % n_colors,
                                              cdict, 1024)


def plot4_seasons(lons, lats, pannels, variables,
                  plot_style='absolute',
                  suptitle='{season}-{variable}',
                  titles=None, units=None,
                  plot_group=None,
                  absolute_kwargs={},
                  anom_kwargs={'cmap': 'RdBu_r', 'cbar_extend': 'both'},
                  pdiff_kwargs={'cmap': 'RdBu_r', 'cbar_extend': 'both',
                                'vmin': -100.0, 'vmax': 100.0},
                  mask=None,
                  figkwds={'figsize': (10, 9)},
                  outdir=None, savekwds={'dpi': 150, 'format': 'png',
                                         'bbox_inches': 'tight'},
                  suptitle_kwargs={'fontsize': 14, 'fontweight': 'bold',
                                   'y': 1.03}):
    """
    Make 4 pannel plots spatial plots.

    There are 5 types of plots availabe:
     - absolute: plot the values in each pannel
     - anomaly: plot the difference of pannels 1, 2, and 3 to pannel 4.
     - percent_difference: plot the percent difference of pannels 1, 2, and 3
        to pannel 4.
     - percent_error: plot the percent error of pannels 1, 2, and 3
        to pannel 4.
     - percent_change: plot the percent change of pannels 1, 2, and 3
        to pannel 4.

    Parameters
    ----------
    pannels : dict
        Dictionary of data pannels with a structure of
        {pannel_number: xray.Dataset} to be accessed as
        [pannel_number][variable][season, y, x].
    variables: str or list like
        Name of variable(s) to plot.
    plot_style : str or list like
        Define which plot styles should be plotted.  Valid values are
        {'absolute', 'anomaly', 'percent_difference', 'percent_error',
        'percent_change'}.  Default is 'absolute'.
    suptitle : str, optional
        String to use as suptitle, the keywords `{season}` and `{variable`}
        will be replaced with the appropriate values using string formating.
    titles : list like, optional
        Names for each pannel, must have length of 4.  The keys of each pannel
            are used by default.
    units : str, optional
        Colorbar label for absolute style plots.  Default is None.
    plot_group : str, optional
        String to use as identifier in output files.  Only used when `outdir`
        is provided.
    absolute_kwargs : dict, optional
        Keyword arguments to pass to `sub_plot_pcolor` when making `absolute`
        sytle plots.
    anom_kwargs : dict, optional
        Keyword arguments to pass to `sub_plot_pcolor` when making `anomaly`
        sytle plots.
    pdiff_kwargs : dict, optional
        Keyword arguments to pass to `sub_plot_pcolor` when making `percent_*`
        sytle plots.
    mask : numpy.ndarray, optional
        2d mask where 0 indicates grid cells that should not be plotted.
    figkwds : dict, optional
        Keyword arguments to pass to `pyplot.figure()`
    outdir : str, optional
        Output directory to write plot files to.  If None, files will not be
        written to files.  Default is None.
    savekwds : dict, optional
        Keyword arguments passed to `pyplot.savefig()`.
    suptitle_kwargs : dict, optional
        Keyword argument passed to `pyplot.suptitle()`.
    """

    if titles is None:
        titles = pannels.keys()
    assert len(titles) == 4

    if 'format' not in savekwds:
        savekwds['format'] = 'png'

    for variable in np.asarray([variables]).squeeze():
        for i, season in enumerate(seasons):
            suptitle = suptitle.format(season=season, variable=variable)

            data = {}
            for j in range(1, 5):
                if mask is None:
                    data[j] = pannels[j][variable].sel(season=season).values
                else:
                    data[j] = np.ma.masked_where(mask == 0,
                                                 pannels[j][variable].sel(
                                                     season=season).values)

            # Calculate the differences in the data pannels
            if 'absolute' in plot_style:
                absolute = True

            if 'anomaly' in plot_style:
                anomaly = {}
                for j in range(1, 4):
                    anomaly[j] = data[j] - data[4]
            else:
                anomaly = False

            if 'percent_difference' in plot_style:
                percent_difference = {}
                for j in range(1, 4):
                    percent_difference[j] = (data[j] - data[4]) / ((data[j] +
                                                                    data[4]) /
                                                                   2.)
            else:
                percent_difference = False

            if 'percent_error' in plot_style:
                percent_error = {}
                for j in range(1, 4):
                    percent_error[j] = (data[j] - data[4]) / data[4]
            else:
                percent_error = False

            if 'percent_change' in plot_style:
                percent_change = {}
                for j in range(1, 4):
                    percent_change[j] = (data[j] - data[4]) / data[j]
            else:
                percent_change = False

            # Make Plots
            if absolute:
                fig = plt.figure(**figkwds)
                for j in range(1, 5):
                    fig.add_subplot(2, 2, j)
                    sub_plot_pcolor(data[j],
                                    title=titles[j], units=units,
                                    **absolute_kwargs)
                fig.suptitle(suptitle, **suptitle_kwargs)
                fig.tight_layout()

                if outdir is not None:
                    outfile = os.path.join(
                        outdir, '%s-%s-%s.%s' % (variable, season, plot_group,
                                                 savekwds['format']))
                    plt.savefig(outfile, **savekwds)

            if anomaly:
                fig = plt.figure(**figkwds)
                for j in range(1, 4):
                    subtitle = '{0}-{1}'.format(titles[j], titles[4])

                    fig.add_subplot(2, 2, j)
                    sub_plot_pcolor(anomaly[j],
                                    title=subtitle, units=units, **anom_kwargs)

                fig.add_subplot(2, 2, 4)
                sub_plot_pcolor(data[4],
                                title=titles[4], units=units,
                                **absolute_kwargs)
                fig.suptitle(suptitle, **suptitle_kwargs)
                fig.tight_layout()

                if outdir is not None:
                    outfile = os.path.join(
                        outdir, '%s-anom-%s-%s.%s' % (variable, season,
                                                      plot_group,
                                                      savekwds['format']))
                    plt.savefig(outfile, **savekwds)

            if percent_difference:
                fig = plt.figure(**figkwds)
                for j in range(1, 4):
                    subtitle = '({0}-{1})/(({0}+{1})/2)'.format(titles[j],
                                                                titles[4])

                    fig.add_subplot(2, 2, j)
                    sub_plot_pcolor(percent_difference[j],
                                    title=subtitle, units='%', **pdiff_kwargs)

                fig.add_subplot(2, 2, 4)
                sub_plot_pcolor(data[4],
                                title=titles[j], units=units,
                                **absolute_kwargs)
                fig.suptitle(suptitle, **suptitle_kwargs)
                fig.tight_layout()

                if outdir is not None:
                    outfile = os.path.join(
                        outdir,
                        '%s-pdiff-%s-%s.%s' % (variable, season, plot_group,
                                               savekwds['format']))
                    plt.savefig(outfile, **savekwds)

            if percent_error:
                fig = plt.figure(**figkwds)
                for j in range(1, 4):
                    subtitle = '({0}-{1})/(({0}+{1})/2)'.format(titles[j],
                                                                titles[4])

                    fig.add_subplot(2, 2, j)
                    sub_plot_pcolor(percent_error[j],
                                    title=subtitle, units='%', **pdiff_kwargs)

                fig.add_subplot(2, 2, 4)
                sub_plot_pcolor(data[4],
                                title=titles[j], units=units,
                                **absolute_kwargs)
                fig.suptitle(suptitle, **suptitle_kwargs)
                fig.tight_layout()

                if outdir is not None:
                    outfile = os.path.join(
                        outdir, '%s-perr-%s-%s.%s' % (variable, season,
                                                      plot_group,
                                                      savekwds['format']))
                    plt.savefig(outfile, **savekwds)

            if percent_change:
                fig = plt.figure(**figkwds)
                for j in range(1, 4):
                    subtitle = '({0}-{1})/(({0}+{1})/2)'.format(titles[j],
                                                                titles[4])

                    fig.add_subplot(2, 2, j)
                    sub_plot_pcolor(percent_change[j],
                                    title=subtitle, units='%', **pdiff_kwargs)

                fig.add_subplot(2, 2, 4)
                sub_plot_pcolor(data[4],
                                title=titles[j], units=units,
                                **absolute_kwargs)
                fig.suptitle(suptitle, **suptitle_kwargs)
                fig.tight_layout()

                if outdir is not None:
                    outfile = os.path.join(
                        outdir, '%s-pchg-%s-%s.%s' % (variable, season,
                                                      plot_group,
                                                      savekwds['format']))
                    plt.savefig(outfile, **savekwds)

default_map = make_bmap()
