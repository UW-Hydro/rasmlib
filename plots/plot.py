#!/usr/bin/env python

""" Plots RASM outputs from netcdf files

The script does the following tasks
 * Reads a configuration file with plot information/settins/paths
 * Fills required settins with defaults if not given in configuration file
 * Makes plots (1, 2, or 4 pannel depending on how many are specified)
 * Saves file to output directory (any existing files will be overwritten)

 An Example of a configuration file is listed below

[Main]
plots: Albedo
outpath: /usr1/jhamman/Dropbox/UW_RASM_shared/r35RB1a_figures/
outformat: png
dpi: 200

[Map_Options]
plot_parallels:-80, 81, 20
plot_meridians:-180, 181, 20

[Projection_Parameters]
projection: npstere
boundinglat: 35
lon_0: -114
lat_ts: 80.5

[Temperature]
plot_title: Surface Temperature
plot_subtitle:Monthly
plot_iterator:01,02,03,04,05,06,07,08,09,10,11,12
plot_anom: True
plot_pannels: 4
plot_cmap: RdBu_r
plot_bounds: False
plot_cbar_extend: both

plot_anom: True
anom_cmap: Jet
anom_bounds: False
anom_cbar_extend: both

pannel0_path: /nfs/hydro6/raid/nijssen/rasm/r33RBVIC60/lnd/hist/r33RBVIC60.vic.ha.1990-1999-[01-12]*.monthly.mean.nc
pannel0_label: r33RBVIC60
pannel0_units: False
pannel0_var: Surft
pannel0_lon: longitude
pannel0_lat: latitude
pannel0_offset: False
pannel0_mult: False
pannel0_div: False
"""

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.cm as GMT
from pylab import axes
import os
from ConfigParser import SafeConfigParser
import argparse

# Default configurations
main_defaults = {'outpath': "./",
                 'outformat': 'png',
                 'dpi': 200}
map_defaults = {'plot_parallels': [-80, 81, 20],
                'plot_meridians': [-180, 181, 20]}
# RASM
projection_defaults = {'projection': 'lcc',
                       'urcrnrlon': 16.90845094934209,
                       'llcrnrlat': 16.534986367884521,
                       'urcrnrlat': 27.511827255753555,
                       'rsphere': 6371200.0,
                       'llcrnrlon': 189.2229322311162,
                       'lon_0': -114,
                       'lat_0': 90}
plot_defaults = {'plot_cmap': 'cm.jet',
                 'plot_cbar_extend': 'both',
                 'plot_anom': True,
                 'anom_cmap': 'cm.RdBu',
                 'anom_cbar_extend': 'both',
                 'anom_type': 'absolute',
                 'plot_bounds': False}


def main():
    Main, Map_Options, Projection_Parameters, Plot_Dict, options = process_command_line()

    # Make sure the outpath exists
    if not os.path.exists(Main['outpath']):
        os.makedirs(Main['outpath'])

    #Loop over all plots in configuration file
    for var in Main['plots']:
        # get list of files from configuration file, loop over
        files,nfiles = get_file_lists(Plot_Dict,var)
        for i in xrange(nfiles):
            Mdata = {}
            Adata = {}
            Coords = {}
            # Get data and attributes for each file and make plots
            for j in xrange(int(Plot_Dict[var]['plot_pannels'])):
                pannel = 'pannel'+str(j)
                if options['verbose']:
                    print 'reading data from', files[pannel][i]
                f = Dataset(files[pannel][i])

                # get offset/multiplier from dictionary
                if Plot_Dict[var][pannel+'_offset']:
                    offset = Plot_Dict[var][pannel+'_offset']
                else:
                    offset = 0.
                if Plot_Dict[var][pannel+'_mult']:
                    mult = Plot_Dict[var][pannel+'_mult']
                else:
                    mult = 1.
                if Plot_Dict[var][pannel+'_div']:
                    div = Plot_Dict[var][pannel+'_div']
                else:
                    div = 1.

                # Get data from netcdf
                var_name = Plot_Dict[var][pannel+'_var']
                lon_name = Plot_Dict[var][pannel+'_lon']
                lat_name = Plot_Dict[var][pannel+'_lat']
                Mdata[pannel+'_var'] = np.squeeze(f.variables[var_name][:])*mult/div + offset
                Coords[pannel+'_lon'] = f.variables[lon_name][:]
                Coords[pannel+'_lat'] = f.variables[lat_name][:]

                # get mask from first dataset and apply to rest.
                if j == 0:
                    mask = Mdata[pannel+'_var'].mask
                else:
                    Mdata[pannel+'_var'] = np.ma.masked_array(Mdata[pannel+'_var'],mask = mask)

                # Check to make sure we have all the necessary attributes
                if not Plot_Dict[var][pannel+'_units']:
                    try:
                        Plot_Dict[var][pannel+'_units'] = f.variables[var_name].units
                    except:
                        Plot_Dict[var][pannel+'_units'] = None
                if not Plot_Dict[var][pannel+'_label']:
                    try:
                        Plot_Dict[var][pannel+'_label'] = f.title
                    except:
                        Plot_Dict[var][pannel+'_label'] = pannel
                f.close()

            # Check to make sure we have a colorbar and title
            if not Plot_Dict[var]['plot_title']:
                Plot_Dict[var]['plot_title'] = var
            if not Plot_Dict[var]['plot_title']:
                Plot_Dict[var]['plot_cmap'] = 'cm.jet'

            # Get colorbar bounds from config or find from data
            bounds = {}
            Mbounds = {}
            if not Plot_Dict[var]['plot_bounds']:
                Mbounds = {}
                Mbounds['vmax'] = max([np.max(Mdata[k]) for k in Mdata])
                Mbounds['vmin'] = min([np.min(Mdata[k]) for k in Mdata])
            else:
                Mbounds['vmin'] = min(Plot_Dict[var]['plot_bounds'])
                Mbounds['vmax'] = max(Plot_Dict[var]['plot_bounds'])
            bounds['Mbounds']= Mbounds

            if options['dryrun']:
                print 'dryrun was sucessful!'
                print 'Plot_Dicts:'
                for var in Plot_Dict:
                    print var,Plot_Dict[var]
                print 'Projection_Parameters:', Projection_Parameters
                print 'Main:',Main
                print 'Map_Options:',Map_Options
                print 'bounds:', bounds
                return

            # Make a plot
            # if Plot_Dict[var]['plot_pannels'] == 1:
            #     plot1()
            # elif Plot_Dict[var]['plot_pannels'] == 2:
            #     plot2()
            # elif Plot_Dict[var]['plot_pannels'] == 3:
            #     plot3()
            if Plot_Dict[var]['plot_pannels'] == 4:
                subtitle = Plot_Dict[var]['plot_iterator'][i]+'-'+Plot_Dict[var]['plot_subtitle']+'-Mean'
                outname = var+'-'+subtitle+'.'+Main['outformat']
                figout = plot4(Mdata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname)
                if options['verbose']:
                    print 'done with ', figout
            else:
                raise ValueError('only have one plot sytle right now')

            # If anomalies are to be plotted, find those now
            if Plot_Dict[var]['plot_anom']:
                Adata = Mdata
                pannels = int(Plot_Dict[var]['plot_pannels'])-1
                for j in xrange(int(pannels)):
                    pannel = 'pannel'+str(j)
                    last = 'pannel'+str(pannels)+'_var'
                    if Plot_Dict[var]['anom_type']=='absolute':
                        Adata[pannel+'_var'] -= Mdata[last]
                    else:
                        Adata[pannel+'_var'] = (Adata[pannel+'_var']-Mdata[last])/Mdata[last]*100

                Abounds = {}
                if not Plot_Dict[var]['anom_bounds']:
                    bmax = max([ma.max(Adata[k]) for k in Adata])
                    bmin = min([ma.min(Adata[k]) for k in Adata])
                    Abounds['vmax'] = np.max(np.absolute([bmax,bmin]))*0.35
                    Abounds['vmin'] = -1*Abounds['vmax']
                else:
                    Abounds['vmin'] = min(Plot_Dict[var]['anom_bounds'])
                    Abounds['vmax'] = max(Plot_Dict[var]['anom_bounds'])
                bounds['Abounds']= Abounds
                # Make a plot
                # if Plot_Dict[var]['plot_pannels'] == 1:
                #     figout = plot1(Adata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname,anom=True)
                # elif Plot_Dict[var]['plot_pannels'] == 2:
                #     figout = plot2(Adata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname,anom=True)
                # elif Plot_Dict[var]['plot_pannels'] == 3:
                #     figout = plot3(Adata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname,anom=True)
                if Plot_Dict[var]['plot_pannels'] == 4:
                    subtitle = Plot_Dict[var]['plot_iterator'][i]+'-'+Plot_Dict[var]['plot_subtitle']+'-Anomaly'
                    outname = var+'-'+subtitle+'.'+Main['outformat']
                    figout = plot4(Adata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname,anom=True)
                    if options['verbose']:
                        print 'done with ', figout
                else:
                    raise ValueError('only have one plot sytle right now')

    return

def process_command_line():
    """
    Read configuration file, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Input Configuration File")
    parser.add_argument("--variable", type=str, help="Limit the plots to a single variable")
    parser.add_argument("--verbose", help="Make script verbose",action='store_true')
    parser.add_argument("--dryrun", help="Do a dryrun without making any plots, check to make sure config file is readable and data exists",action='store_true')
    args = parser.parse_args()

    options = {}
    options['verbose'] = args.verbose
    options['dryrun'] = args.dryrun

    config_dict = read_config(args.config)

    # unwrap the config dict, and check to make sure all the default keys are present
    Main = config_dict.pop('Main')
    for key in main_defaults:
        if key not in Main:
            Main[key] = main_defaults[key]
    if not isinstance(Main['plots'], list):
        Main['plots'] = list([Main['plots']])
    if args.variable:
        Main['plots'] = [args.variable]


    Map_Options = config_dict.pop('Map_Options')
    for key in map_defaults:
        if key not in Map_Options:
            Map_Options[key] = map_defaults[key]

    Projection_Parameters = config_dict.pop('Projection_Parameters')
    for key in projection_defaults:
        if key not in Projection_Parameters:
            Projection_Parameters[key] = projection_defaults[key]

    # Deal with plots section
    for var in Main['plots']:
        print config_dict[var]['plot_pannels']
        for num in xrange(int(config_dict[var]['plot_pannels'])):
            pannel = 'pannel'+str(num)
            if pannel+'_path' not in config_dict[var]:
                raise IOError('Need input file for '+pannel)
            if pannel+'_label' not in config_dict[var]:
                config_dict[var][pannel+'_label'] = False
            if pannel+'_units' not in config_dict[var]:
                config_dict[var][pannel+'_units'] = False
            if pannel+'_var' not in config_dict[var][pannel+'_var']:
                config_dict[var][pannel+'_var'] = var
            if pannel+'_lon' not in config_dict[var]:
                config_dict[var][pannel+'_lon'] = 'longitude'
            if pannel+'_lat' not in config_dict[var]:
                config_dict[var][pannel+'_lat'] = 'latitude'
            if pannel+'_offset' not in config_dict[var]:
                config_dict[var][pannel+'_offset'] = False
            if pannel+'_mult' not in config_dict[var]:
                config_dict[var][pannel+'_mult'] = False
            if pannel+'_div' not in config_dict[var]:
                config_dict[var][pannel+'_div'] = False

    return Main, Map_Options, Projection_Parameters, config_dict, options

def get_file_lists(Plot_Dict,var):
    """
    Get a list of files based on the wildcard (if any) from the configuration file.
    """
    files = {}
    for i in xrange(int(Plot_Dict[var]['plot_pannels'])):
        pannel = 'pannel'+str(i)
        files[pannel] = []
        for j in Plot_Dict[var]['plot_iterator']:
            files[pannel].append(Plot_Dict[var][pannel+'_path'].replace('***',j))
        nfiles = len(files[pannel])

    return files,nfiles

# def plot1(data,coords,bounds,Main,Map_Options,Projection_Parameters,Var_Dict):
#     """
#     Make a single pannel plot
#     """
#     fig = plt.figure()
#     pannel = 'pannel0'
#     units = Var_Dict[pannel+'_units']
#     m = Basemap(**Projection_Parameters)
#     xi,yi = m(coords[pannel+'_lon'], coords[pannel+'_lat'])
#     m.drawcoastlines()
#     m.drawparallels(np.arange(*Map_Options['plot_parallels']))
#     m.drawmeridians(np.arange(*Map_Options['plot_meridians']))
#     cs = m.pcolor(xi,yi, data[pannel+'_var'])

#     plt.show()

#     return

# def plot2in4():
#     """
#     Make a 2 variables in a 4 pannel plot
#     """


#     return

# def plot3(data,coords,bounds,Main,Map_Options,Projection_Parameters,Var_Dict,subtitle,outname,anom=False):
#     """
#     Make a 3 pannel  horizontal plot
#     """
#     fig = plt.figure()
#     b = 0.04; c = 0;
#     b = b/3; c = 0.5*c;
#     xxs = [0.15+b, 0.33-b, 0.67-b]
#     yys = [0.125+c, 0.125+c, 0.125+c]
#     subwidth = 0.35; subheight = 0.33

#     for i in xrange(3):
#         pannel = 'pannel'+str(i)
#         units = Var_Dict[pannel+'_units']
#         label = Var_Dict[pannel+'_label']
#         axes([xxs[i], yys[i], subwidth, subheight])

#         m = Basemap(**Projection_Parameters)
#         m.drawcoastlines()
#         m.drawparallels(np.arange(*Map_Options['plot_parallels']))
#         m.drawmeridians(np.arange(*Map_Options['plot_meridians']))

#         xi,yi = m(coords[pannel+'_lon'], coords[pannel+'_lat'])

#         if anom and i==3:
#             cmap =  cmap_discretize(Var_Dict['plot_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Mbounds'])
#             extend = Var_Dict['plot_cbar_extend']
#             title = Var_Dict[pannel+'_label']

#         elif anom:
#             cmap =  cmap_discretize(Var_Dict['anom_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Abounds'])
#             extend = Var_Dict['anom_cbar_extend']
#             title = Var_Dict[pannel+'_label']+' - '+Var_Dict['pannel3_label']
#         else:
#             cmap =  cmap_discretize(Var_Dict['plot_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Mbounds'])
#             extend = Var_Dict['plot_cbar_extend']
#             title = Var_Dict[pannel+'_label']

#         if i%2 == 0:
#             cbar = m.colorbar(cs, location='left',extend = extend )
#             cbar.ax.yaxis.set_label_position('left')
#             cbar.ax.yaxis.set_ticks_position('left')
#         else:
#             cbar = m.colorbar(cs,location='right',extend=extend)
#         cbar.set_label(units)
#         for t in cbar.ax.get_yticklabels():
#             t.set_fontsize(9)

#         plt.title(title, fontsize=12)

#     fig.suptitle(Var_Dict['plot_title'], fontsize=14, fontweight='bold')
#     fig.text(0.5, 0.92, subtitle, ha='center',fontsize=12)
#     figout = os.path.join(Main['outpath'],outname)
#     plt.savefig(figout, format=Main['outformat'], dpi=Main['dpi'],bbox_inches='tight', pad_inches=0.1)
#     plt.close()
#     return figout

# def plot4(data,coords,bounds,Main,Map_Options,Projection_Parameters,Var_Dict,subtitle,outname,anom=False):
#     """
#     Make a 4 pannel plot
#     """
#     fig = plt.figure()
#     b = 0.04; c = 0;
#     b = 0.5*b;c = 0.5*c;
#     xxs = [0.15+b, 0.5-b, 0.15+b, 0.5-b]
#     yys = [0.525-c, 0.525-c, 0.125+c, 0.125+c]
#     subwidth = 0.35; subheight = 0.33

#     for i in xrange(4):
#         pannel = 'pannel'+str(i)
#         units = Var_Dict[pannel+'_units']
#         label = Var_Dict[pannel+'_label']
#         axes([xxs[i], yys[i], subwidth, subheight])

#         m = Basemap(**Projection_Parameters)
#         m.drawcoastlines()
#         m.drawparallels(np.arange(*Map_Options['plot_parallels']))
#         m.drawmeridians(np.arange(*Map_Options['plot_meridians']))

#         xi,yi = m(coords[pannel+'_lon'], coords[pannel+'_lat'])

#         if anom and i==3:
#             cmap =  cmap_discretize(Var_Dict['plot_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Mbounds'])
#             extend = Var_Dict['plot_cbar_extend']
#             title = Var_Dict[pannel+'_label']

#         elif anom:
#             cmap =  cmap_discretize(Var_Dict['anom_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Abounds'])
#             extend = Var_Dict['anom_cbar_extend']
#             title = Var_Dict[pannel+'_label']+' - '+Var_Dict['pannel3_label']
#         else:
#             cmap =  cmap_discretize(Var_Dict['plot_cmap'])
#             cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cmap,**bounds['Mbounds'])
#             extend = Var_Dict['plot_cbar_extend']
#             title = Var_Dict[pannel+'_label']

#         if i%2 == 0:
#             cbar = m.colorbar(cs, location='left',extend = extend )
#             cbar.ax.yaxis.set_label_position('left')
#             cbar.ax.yaxis.set_ticks_position('left')
#         else:
#             cbar = m.colorbar(cs,location='right',extend=extend)
#         cbar.set_label(units)
#         for t in cbar.ax.get_yticklabels():
#             t.set_fontsize(9)

#         plt.title(title, fontsize=12)

#     fig.suptitle(Var_Dict['plot_title'], fontsize=14, fontweight='bold')
#     fig.text(0.5, 0.92, subtitle, ha='center',fontsize=12)
#     figout = os.path.join(Main['outpath'],outname)
#     plt.savefig(figout, format=Main['outformat'], dpi=Main['dpi'],bbox_inches='tight', pad_inches=0.1)
#     plt.close()
#     return figout

def setup_map(lons, lats, **kwargs):

    proj_params = {'projection': 'lcc',
                   'llcrnrlat': lats[0, 0],
                   'urcrnrlat': lats[-1, -1],
                   'llcrnrlon': lons[0, 0],
                   'urcrnrlon': lons[-1, -1],
                   'rsphere': 6371200.0,
                   'lon_0': -114,
                   'lat_0': 90}

    m = Basemap(**proj_params)
    xi,yi = m(np.squeeze(lons),np.squeeze(lats))

    return m, xi, yi

def sub_plot_pcolor(m, xi, yi, data, title='', cmap=None,
                    cbar=True, cbar_location='bottom', cbar_extend=False,
                    units=None, vmin=None, vmax=None):
    """individual subplot function"""

    sp = m.pcolormesh(xi,yi,np.squeeze(data), vmin=vmin, vmax=vmax, cmap=cmap)

    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawcoastlines(color='k',linewidth=0.25);

    if title:
        plt.title(title, size=14)
    if cbar:
        cbar = m.colorbar(location=cbar_location, extend=cbar_extend, cmap=cmap)
    if (units and cbar):
        cbar.set_label(units)

    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(9)

    return sp



def plot4(data, coords, bounds, Main, Map_Options, Projection_Parameters, Var_Dict, subtitle, outname, anom=False):
    """
    Make a 4 pannel plot
    """
    fig = plt.figure(figsize=(10,9.5))

    for i in xrange(4):
        ax = fig.add_subplot(2, 2, i)
        pannel = 'pannel'+str(i)
        units = Var_Dict[pannel+'_units']
        label = Var_Dict[pannel+'_label']

        m, xi, yi = setup_map(coords[pannel+'_lon'], coords[pannel+'_lat'], **Projection_Parameters)

        m.drawcoastlines()
        m.drawparallels(np.arange(*Map_Options['plot_parallels']))
        m.drawmeridians(np.arange(*Map_Options['plot_meridians']))

        if not anom or i==3:
            cmap =  cmap_discretize(Var_Dict['plot_cmap'])
            cs = sub_plot_pcolor(m, xi, yi, data[pannel+'_var'],
                                 title=label, cmap=cmap, units=units,
                                 cbar_extend=Var_Dict['plot_cbar_extend'],
                                 vmin=bounds['Mbounds']['vmin'],
                                 vmax=bounds['Mbounds']['vmax'])
        else:
            cmap =  cmap_discretize(Var_Dict['anom_cmap'])
            if Var_Dict['anom_type']=='absolute':
                title = label+' - '+Var_Dict['pannel3_label']
            else:
                title = "%%diff: (%s - %s)/ %s" %(label, Var_Dict['pannel3_label'], Var_Dict['pannel3_label'])
            cs = sub_plot_pcolor(m, xi, yi, data[pannel+'_var'],
                                 title=title, cmap=cmap, units=units,
                                 cbar_extend=Var_Dict['anom_cbar_extend'],
                                 vmin=bounds['Abounds']['vmin'],
                                 vmax=bounds['Abounds']['vmax'])

    fig.suptitle(Var_Dict['plot_title'], fontsize=14, fontweight='bold', y=1.05)
    fig.text(0.5, 1, subtitle, ha='center', fontsize=12)
    fig.tight_layout()
    figout = os.path.join(Main['outpath'], outname)

    plt.savefig(figout, format=Main['outformat'], dpi=Main['dpi'], bbox_inches='tight', pad_inches=0.1)
    plt.close()
    return figout

def cmap_discretize(cmap, N=10):

     cmap = cm.get_cmap(eval(cmap))
     colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
     colors_rgba = cmap(colors_i)
     indices = np.linspace(0, 1., N+1)
     cdict = {}
     for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1)]

     return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

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

if __name__ == "__main__":
    main()
