"""
rasmlib plot utilities
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np
from ..calendar import seasons
from .climatology import season_mean, annual_mean
from scipy.stats import ttest_ind

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


def cmap_discretize_2(cmap, n_colors):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        n_colors: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n_colors + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n_colors + 1)]
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d" % n_colors,
                                              cdict, 1024)


def colorbar_index(ncolors, cmap, ticklabels=None, cbar_shrink=1.0):
    cmap = cmap_discretize_2(cmap, ncolors)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors + 0.5)
    colorbar = plt.colorbar(mappable, shrink=cbar_shrink, pad=0.01)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    if ticklabels is None:
        colorbar.set_ticklabels(range(ncolors))
    else:
        colorbar.set_ticklabels(ticklabels)
    return colorbar


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
    return


def plot_zonal_means(lats, data, bounds=None, mask=None, weights=None,
                     num_bins=100, xlabel='Latitude', ylabel='', title=''):
    """ """
    sns.set_style("darkgrid", {"grid.linewidth": .5, "axes.facecolor": ".9"})

    plt.figure()

    if mask is None:
        mask = np.ones_like(lats)
    if bounds is None:
        bounds = [lats.min(), lats.max()]
    bins = np.linspace(*bounds, num=num_bins + 1)
    lat_means = np.zeros(num_bins)
    data_means = {}
    for k, v in data.items():
        data_means[k] = np.zeros(num_bins)

    for i in range(num_bins):
        y0 = bins[i]
        y1 = bins[i + 1]
        lat_means[i] = bins[i:i + 1].mean()
        ys, xs = np.nonzero((mask != 0) & (lats >= y0) & (lats < y1))
        if weights is None:
            for k, d in data.items():
                data_means[k][i] = np.nanmean(d[ys, xs])
        else:
            for k, d in data.items():
                data_means[k][i] = (np.nansum(d[ys, xs] * weights[ys, xs]) /
                                    np.nansum(weights[ys, xs]))

    for k, d in data_means.items():
        plt.plot(lat_means, d, label=k)

    plt.legend()
    plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
    return


def make_plot_data(ds1, ds2, ny, nx, mask=None, relative=False,
                   start=None, end=None):

    ncols = 5
    nrows = 3

    plot_data = np.ma.empty((nrows, ncols, ny, nx))
    hatch_data = np.ma.empty((ncols, ny, nx))

    # ds1
    temp1 = season_mean(ds1.sel(time=slice(start, end)),
                        calendar='noleap').squeeze()
    plot_data[0, 4] = annual_mean(ds1.sel(time=slice(start, end)),
                                  calendar='noleap').squeeze().values

    # ds2
    temp2 = season_mean(ds2.sel(time=slice(start, end)),
                        calendar='standard').squeeze()
    plot_data[1, 4] = annual_mean(ds2.sel(time=slice(start, end)),
                                  calendar='standard').squeeze().values

    for i, season in enumerate(seasons):
        plot_data[0, i] = temp1.sel(season=season).values
        plot_data[1, i] = temp2.sel(season=season).values

    # diff
    if relative:
        # percent difference
        plot_data[2] = (100. * (plot_data[0] - plot_data[1]) /
                        ((plot_data[0] + plot_data[1]) / 2))
    else:
        plot_data[2] = plot_data[0] - plot_data[1]

    # hatch_data
    ds1groups = dict(ds1.groupby('time.season'))
    ds2groups = dict(ds2.groupby('time.season'))
    for i, season in enumerate(seasons):
        _, hatch_data[i] = ttest_ind(ds1groups[season].values.squeeze(),
                                     ds2groups[season].values.squeeze(),
                                     equal_var=False)
    _, hatch_data[4] = ttest_ind(ds1.values.squeeze(),
                                 ds2.values.squeeze(),
                                 equal_var=False)

    # mask
    if mask is not None:
        temp = np.ma.masked_where(mask == 0, mask)
        mask = np.broadcast_arrays(plot_data, temp.mask)[1]
        plot_data = np.ma.masked_where(mask, plot_data)

    return plot_data, hatch_data


def plot2(plot_data,
          cmap='Spectral_r',
          amap='RdBu_r',
          ylabels=None,
          titles=('DJF', 'MAM', 'JJA', 'SON', 'ANNUAL'),
          suptitle=None,
          cbar_label='',
          abar_label='',
          vmin=None,
          vmax=None,
          amin=None,
          amax=None,
          cbar_extend='neither',
          abar_extend='neither',
          stack_cbars=False,
          map_obj=default_map,
          anom_hatching=None):
    """"""
    ncols = 5
    nrows = 3

    # Make sure plot_data is the right size
    assert(map_obj.xi.shape == map_obj.yi.shape)
    assert(map_obj.xi.shape == plot_data.shape[2:])
    assert(plot_data.shape[:2] == (nrows, ncols))

    # set data ranges
    if vmin is None:
        vmin = plot_data[:2].min()
    if vmax is None:
        vmax = plot_data[:2].max()
    if amin is None:
        amin = plot_data[2].min()
    if amax is None:
        amax = plot_data[2].max()

    # Set colorbar norms and ticks
    assert(type(cmap) == str)
    cn = 10
    cmap = cmap_discretize(cmap, n_colors=cn)
    cnorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cticks = np.linspace(vmin, vmax, num=cn + 1)
    assert(type(amap) == str)
    an = 10
    amap = cmap_discretize(amap, n_colors=an)
    anorm = mpl.colors.Normalize(vmin=amin, vmax=amax)
    aticks = np.linspace(amin, amax, num=an + 1)

    if abar_label and not cbar_label:
        cbar_label = abar_label

    # Copy colormap and data ranges to iterables
    cmaps = (cmap, cmap, amap)
    vmins = (vmin, vmin, amin)
    vmaxs = (vmax, vmax, amax)

    # Make the plot
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(11, 5.5))
    for (i, j), ax in np.ndenumerate(axes):
        plt.sca(ax)
        sub_plot_pcolor(plot_data[i, j],
                        cmap=cmaps[i],
                        cbar=None,
                        vmin=vmins[i],
                        vmax=vmaxs[i],
                        map_obj=map_obj,
                        ax=ax)
        if i == 0 and titles is not None:
            ax.set_title(titles[j])
        if j == 0 and ylabels is not None:
            ax.set_ylabel(ylabels[i])
        if i == nrows - 1 and anom_hatching is not None:
            map_obj.m.contourf(map_obj.xi, map_obj.yi,
                               100 * (anom_hatching[j]),
                               [0, 5],
                               cmap=plt.get_cmap('gray'),
                               hatches=['....', None],
                               alpha=0,
                               ax=ax)

    # Add figure title
    if suptitle is not None:
        fig.suptitle(suptitle, fontsize=16, fontweight='roman', y=1.02)
    plt.tight_layout()

    # Color bars
    if stack_cbars:
        ax1 = fig.add_axes([0.995, 0.37, 0.015, 0.54])
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=cnorm,
                                        orientation='vertical',
                                        extend=cbar_extend,
                                        ticks=cticks)
        ax2 = fig.add_axes([0.995, 0.08, 0.015, 0.23])
        cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=amap, norm=anorm,
                                        orientation='vertical',
                                        extend=abar_extend,
                                        ticks=aticks,
                                        extendfrac=0.12)
        if cbar_label:
            cb1.set_label(cbar_label, rotation=90)
            cb2.set_label(abar_label, rotation=90)

    else:
        ax1 = fig.add_axes([0.995, 0.08, 0.015, 0.83])
        cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=cnorm,
                                        orientation='vertical',
                                        extend=cbar_extend,
                                        ticks=cticks)
        ax2 = fig.add_axes([1.05, 0.08, 0.015, 0.83])
        cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=amap, norm=anorm,
                                        orientation='vertical',
                                        extend=abar_extend,
                                        ticks=aticks)
        if cbar_label:
            cb1.set_label(cbar_label, y=0.005, labelpad=-10, rotation=0)
        if abar_label:
            cb2.set_label(abar_label, y=0.005, labelpad=-10, rotation=0)

    return fig, axes


def plot_n(monthly_means,
           annual_means=False,
           cmap='Spectral_r',
           amap='RdBu_r',
           vmin=None,
           vmax=None,
           amin=None,
           amax=None,
           map_obj=default_map,
           cbar_label='',
           abar_label='',
           cbar_extend='neither',
           abar_extend='neither',
           mask=None):

    ''''''

    nrows = len(monthly_means)
    if annual_means:
        ncols = 5
    else:
        ncols = 4

    width = 11
    height = 1.55 * nrows + 0.6

    # Set colorbar norms and ticks
    assert(type(cmap) == str)
    cn = 10
    cmap = cmap_discretize(cmap, n_colors=cn)
    cnorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cticks = np.linspace(vmin, vmax, num=cn + 1)
    assert(type(amap) == str)
    an = 11
    amap = cmap_discretize(amap, n_colors=an)
    anorm = mpl.colors.Normalize(vmin=amin, vmax=amax)
    aticks = np.linspace(amin, amax, num=an + 1)

    if abar_label and not cbar_label:
        cbar_label = abar_label

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width, height),
                             squeeze=False)

    plt.subplots_adjust(left=0.125, bottom=0.05, right=0.9, top=0.9,
                        wspace=0.05, hspace=0.05)

    for i, (k, ds) in enumerate(monthly_means.items()):

        if i == 0:
            i0_ylabel = k
            axes[i, 0].set_ylabel(i0_ylabel)
            seas_mean0 = season_mean(ds)
            for j, season in enumerate(seasons):
                plt.sca(axes[i, j])
                sub_plot_pcolor(
                    np.ma.masked_where(
                        mask, seas_mean0.sel(season=season).values.squeeze()),
                    cmap=cmap,
                    map_obj=map_obj,
                    ax=axes[i, j],
                    cbar=None,
                    vmin=vmin,
                    vmax=vmax)
                axes[i, j].set_title(season)
            if annual_means:
                ann_mean0 = annual_mean(ds)
                plt.sca(axes[i, 4])
                sub_plot_pcolor(
                    np.ma.masked_where(mask, ann_mean0.values.squeeze()),
                    cmap=cmap,
                    map_obj=map_obj,
                    ax=axes[i, 4],
                    cbar=None,
                    vmin=vmin,
                    vmax=vmax)
                axes[i, 4].set_title('ANNUAL')
        else:
            axes[i, 0].set_ylabel('{0}\n — {1}'.format(i0_ylabel, k))
            seas_meani = seas_mean0 - season_mean(ds)
            for j, season in enumerate(seasons):
                plt.sca(axes[i, j])
                sub_plot_pcolor(
                    np.ma.masked_where(
                        mask, seas_meani.sel(season=season).values.squeeze()),
                    cmap=amap,
                    map_obj=map_obj,
                    ax=axes[i, j],
                    cbar=None,
                    vmin=amin,
                    vmax=amax)
            if annual_means:
                ann_meani = ann_mean0 - annual_mean(ds)
                plt.sca(axes[i, 4])
                sub_plot_pcolor(
                    np.ma.masked_where(mask, ann_meani.values.squeeze()),
                    cmap=amap,
                    map_obj=map_obj,
                    ax=axes[i, 4],
                    cbar=None,
                    vmin=amin,
                    vmax=amax)

    ax1 = fig.add_axes([0.05, 0.0, 0.44, 0.015])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=cnorm,
                                    orientation='horizontal',
                                    extend=cbar_extend,
                                    ticks=cticks)
    ax2 = fig.add_axes([0.55, 0.0, 0.44, 0.015])
    cb2 = mpl.colorbar.ColorbarBase(ax2, cmap=amap, norm=anorm,
                                    orientation='horizontal',
                                    extend=abar_extend,
                                    ticks=aticks)
    if cbar_label:
        cb1.set_label(cbar_label)
    if abar_label:
        cb2.set_label(abar_label)

    plt.tight_layout(pad=.95)

    return fig, axes


def plot_percentile_and_mean(x, data, percentiles=(25, 75), ax=None,
                             fill_kwargs={}, plot_kwargs={}):
    """plot percentile and mean"""
    if ax is None:
        ax = plt.gca()
    p1, p2 = np.nanpercentile(data, percentiles, axis=0)
    ax.fill_between(x, p1, p2, **fill_kwargs)
    ax.plot(x, np.nanmean(data, axis=0), **plot_kwargs)

default_map = make_bmap()
