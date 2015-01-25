import numpy as np
import matplotlib
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
import matplotlib.pyplot as plt
from ConfigParser import SafeConfigParser


def subplot_pcolor(lons, lats, data, title=None, cmap=cm.jet,
                   vmin=None, vmax=None, cbar=True, cbar_location='bottom',
                   units=None, projection=None):

    if not vmin:
        vmin = data.min()
    if not vmax:
        vmax = data.max()

    if not projection:
        projection = define_projection(lons, lats)

    m = Basemap(**projection)
    xi, yi = m(np.squeeze(lons), np.squeeze(lats))
    sp = m.pcolormesh(xi, yi, np.squeeze(data), vmin=vmin, vmax=vmax,
                      cmap=cmap)
    m.drawparallels(np.arange(-80., 81., 20.))
    m.drawmeridians(np.arange(-180., 181., 20.))
    m.drawcoastlines(color='k', linewidth=0.25)

    if title:
        plt.title(title, size=13)
    if cbar:
        cbar = m.colorbar(location=cbar_location)
    cbar.set_label(units)

    return sp


def cmap_discretize(cmap, N=10):

    cmap = cm.get_cmap(eval(cmap))
    colors_i = np.concatenate((np.linspace(0,  1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1, ki], colors_rgba[i, ki])
                      for i in xrange(N+1)]

    return matplotlib.colors.LinearSegmentedColormap('{0}_{1}'.format(cmap.name, N),
                                                     cdict, 1024)


def define_projection(lons, lats):

    # note:  need to add check so that lats/lons are 2d and ordered correctly
    # default to lambert conformal basemap projection
    projection = {'urcrnrlat': lats[0, -1],
                  'urcrnrlon': lons[0, -1],
                  'llcrnrlat': lats[-1, 0],
                  'llcrnrlon': lons[-1, 0],
                  'projection': 'lcc',
                  'rsphere': 6371200.0,
                  'lon_0':  lons.mean(),
                  'lat_0': lats.mean()}
    return projection


# -------------------------------------------------------------------- #
# Read the Configuration File
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
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
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
                return float(value)
            except:
                return value
    else:
        try:
            return map(float, val_list)
        except:
            return val_list
# -------------------------------------------------------------------- #
