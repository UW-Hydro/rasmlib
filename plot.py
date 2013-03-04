#!/usr/local/bin/python

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

[Projection_Parmaeters]
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

"""

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from pylab import axes
from matplotlib import rc
import pdb
import os, glob
import ConfigParser
import argparse

def main():
    Main,Map_Options,Projection_Parameters,Plot_Dict = process_command_line()
    
    #Loop over all plots in configuration file
    for var in Main['plots']:
        # get list of files from configuration file, loop over
        files,nfiles = get_file_lists(Plot_Dict,var)
        for i in xrange(nfiles):
            Mdata = {}
            Adata = {}
            Coords = {}
            # Get data and attributes for each file and make plots
            for j in xrange(Plot_Dict[var]['plot_pannels']):
                pannel = 'pannel'+str(j)
                print 'reading data from', files[pannel][i]
                f = Dataset(files[pannel][i])

                # get offset/multiplier from dictionary
                if Plot_Dict[var][pannel+'_offset']:
                    offset = Plot_Dict[var][pannel+'_offset']
                else:
                    offset = 0
                if Plot_Dict[var][pannel+'_mult']:
                    mult = Plot_Dict[var][pannel+'_mult']
                else:
                    mult = 1

                # Get data from netcdf 
                var_name = Plot_Dict[var][pannel+'_var']
                lon_name = Plot_Dict[var][pannel+'_lon']
                lat_name = Plot_Dict[var][pannel+'_lat']
                Mdata[pannel+'_var'] = np.squeeze(f.variables[var_name][:])*mult + offset
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
                Plot_Dict[var]['plot_cmap'] = 'jet'

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
                

            # Make a plot
            if Plot_Dict[var]['plot_pannels'] == 1:
                plot1()
            elif Plot_Dict[var]['plot_pannels'] == 2:
                plot2()
            elif Plot_Dict[var]['plot_pannels'] == 3:
                plot3()
            elif Plot_Dict[var]['plot_pannels'] == 4:
                subtitle = Plot_Dict[var]['plot_iterator'][i]+'-'+Plot_Dict[var]['plot_subtitle']+'-Mean'
                outname = var+'-'+subtitle+'.'+Main['outformat']
                plot4(Mdata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname)

            # If anomalies are to be plotted, find those now
            if Plot_Dict[var]['plot_anom']:
                Adata = Mdata
                pannels = Plot_Dict[var]['plot_pannels']-1
                for j in xrange(pannels):
                    pannel = 'pannel'+str(j)
                    last = 'pannel'+str(pannels)+'_var'
                    Adata[pannel+'_var'] -= Mdata[last]

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
                if Plot_Dict[var]['plot_pannels'] == 1:
                    plot1()
                elif Plot_Dict[var]['plot_pannels'] == 2:
                    plot2()
                elif Plot_Dict[var]['plot_pannels'] == 3:
                    plot3()
                elif Plot_Dict[var]['plot_pannels'] == 4:
                    subtitle = Plot_Dict[var]['plot_iterator'][i]+'-'+Plot_Dict[var]['plot_subtitle']+'-Anomaly'
                    outname = var+'-'+subtitle+'.'+Main['outformat']
                    plot4(Adata,Coords,bounds,Main,Map_Options,Projection_Parameters,Plot_Dict[var],subtitle,outname,anom=True)

    return
            
            
 
def process_command_line():
    """
    Read configuration file, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="Input Configuration File")
    args = parser.parse_args()
    if args.config:
        configFile = args.config
    else:
        raise IOError('Need to supply configuration file, ./makeplot.py configfile.cfg, type ./makeplot.py -h for more info')
    config = ConfigParser.ConfigParser()
    config.read(configFile)
    config.sections()

    # Parse Main Section
    Main = {}
    Main['plots'] = config.get('Main','plots').split(',')
    Main['outpath'] = config.get('Main','outpath')
    Main['outformat'] = config.get('Main','outformat')
    Main['dpi'] = int(config.get('Main','dpi'))

    # Parse Map Options
    Map_Options = {}
    Map_Options['plot_parallels'] = map(float,config.get('Map_Options','plot_parallels').split(','))
    Map_Options['plot_meridians'] = map(float,config.get('Map_Options','plot_meridians').split(','))

    # Parse Projection Parameters
    Projection_Parameters = {}
    Projection_Parameters['projection'] = config.get('Projection_Parmaeters','projection')
    Projection_Parameters['boundinglat'] = float(config.get('Projection_Parmaeters','boundinglat'))
    Projection_Parameters['lon_0'] = float(config.get('Projection_Parmaeters','lon_0'))
    Projection_Parameters['lat_ts'] = float(config.get('Projection_Parmaeters','lat_ts'))

    # Parse Plots
    Plot_Dict = {}
    for var in Main['plots']:
        Plot_Dict[var] = {}
        Plot_Dict[var]['plot_title'] = config.get(var,'plot_title')
        Plot_Dict[var]['plot_subtitle'] = config.get(var,'plot_subtitle')
        Plot_Dict[var]['plot_iterator'] = config.get(var,'plot_iterator').split(',')
        Plot_Dict[var]['plot_anom'] = config.getboolean(var,'plot_anom')
        Plot_Dict[var]['plot_pannels'] = int(config.get(var,'plot_pannels'))
        try:
            Plot_Dict[var]['plot_bounds'] = map(float,config.get(var,'plot_bounds').split(','))
        except:
            Plot_Dict[var]['plot_bounds'] = False
        try:
            Plot_Dict[var]['plot_cmap'] = config.get(var,'plot_cmap')
        except:
            Plot_Dict[var]['plot_cmap'] = 'RdBu_r'
        try:
            Plot_Dict[var]['plot_cbar_extend'] = config.get(var,'plot_cbar_extend')
        except:
             Plot_Dict[var]['plot_cbar_extend'] = 'neither'
        if  Plot_Dict[var]['plot_anom']:
            try:
                Plot_Dict[var]['anom_cmap'] = config.get(var,'anom_cmap')
            except:
                Plot_Dict[var]['anom_cmap'] = 'Jet'
            try:
                Plot_Dict[var]['anom_bounds'] = map(float,config.get(var,'anom_bounds').split(','))
            except:
                Plot_Dict[var]['anom_bounds'] = False
            try:
                Plot_Dict[var]['anom_cbar_extend'] = config.get(var,'anom_cbar_extend')
            except:
                Plot_Dict[var]['anom_cbar_extend'] = 'neither'
        for num in xrange(int(config.get(var,'plot_pannels'))):
            pannel = 'pannel'+str(num)
            try:
                Plot_Dict[var][pannel+'_path'] = config.get(var,pannel+'_path')
            except:
                raise IOError('Need input file for '+pannel)
            try:
                Plot_Dict[var][pannel+'_label'] = config.get(var,pannel+'_label')
            except:
                Plot_Dict[var][pannel+'_label'] = False
            try:
                Plot_Dict[var][pannel+'_units'] = config.get(var,pannel+'_units') 
            except:
                Plot_Dict[var][pannel+'_units'] = False
            try:
                Plot_Dict[var][pannel+'_var'] = config.get(var,pannel+'_var')
            except:
                Plot_Dict[var][pannel+'_var'] = var
            try:
                Plot_Dict[var][pannel+'_lon'] = config.get(var,pannel+'_lon')
            except:
                Plot_Dict[var][pannel+'_lon'] = 'longitude'
            try:
                Plot_Dict[var][pannel+'_lat'] = config.get(var,pannel+'_lat')
            except:
                Plot_Dict[var][pannel+'_lat'] = 'latitude'
            try:
                Plot_Dict[var][pannel+'_offset'] = float(config.get(var,pannel+'_offset'))
            except:
                Plot_Dict[var][pannel+'_offset'] = False
            try:
                Plot_Dict[var][pannel+'_mult'] = float(config.get(var,pannel+'_mult'))
            except:
                Plot_Dict[var][pannel+'_mult'] = False
    return Main,Map_Options,Projection_Parameters,Plot_Dict

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

def plot1(data,coords,bounds,Main,Map_Options,Projection_Parameters,Var_Dict):
    """
    Make a single pannel plot
    """
    fig = plt.figure()
    pannel = 'pannel0'
    units = Var_Dict[pannel+'_units']
    m = Basemap(**Projection_Parameters)
    xi,yi = m(coords[pannel+'_lon'], coords[pannel+'_lat'])
    m.drawcoastlines()
    m.drawparallels(np.arange(*Map_Options['plot_parallels']))
    m.drawmeridians(np.arange(*Map_Options['plot_meridians']))
    cs = m.pcolor(xi,yi, data[pannel+'_var'])
    
    plt.show()

    return

def plot2():
    """
    Make a 2 pannel plot
    """
    return

def plot3():
    """
    Make a 3 pannel plot
    """
    return

def plot4(data,coords,bounds,Main,Map_Options,Projection_Parameters,Var_Dict,subtitle,outname,anom=False):
    """
    Make a 4 pannel plot
    """
    fig = plt.figure()
    b = 0.04; c = 0;
    b = 0.5*b;c = 0.5*c;
    xxs = [0.15+b, 0.5-b, 0.15+b, 0.5-b]
    yys = [0.525-c, 0.525-c, 0.125+c, 0.125+c]
    subwidth = 0.35; subheight = 0.33
    
    for i in xrange(4):
        pannel = 'pannel'+str(i)
        units = Var_Dict[pannel+'_units']
        label = Var_Dict[pannel+'_label']
        axes([xxs[i], yys[i], subwidth, subheight])
        
        m = Basemap(**Projection_Parameters)
        m.drawcoastlines()
        m.drawparallels(np.arange(*Map_Options['plot_parallels']))
        m.drawmeridians(np.arange(*Map_Options['plot_meridians']))

        xi,yi = m(coords[pannel+'_lon'], coords[pannel+'_lat'])

        if anom and i==3:
            cm = matplotlib.cm.get_cmap(Var_Dict['plot_cmap'])
            cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cm,**bounds['Mbounds'])
            extend = Var_Dict['plot_cbar_extend']
        elif anom:
            cm = matplotlib.cm.get_cmap(Var_Dict['anom_cmap'])
            cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cm,**bounds['Abounds'])
            extend = Var_Dict['anom_cbar_extend']
        else:
            cm = matplotlib.cm.get_cmap(Var_Dict['plot_cmap'])
            cs = m.pcolor(xi,yi, data[pannel+'_var'],cmap=cm,**bounds['Mbounds'])
            extend = Var_Dict['plot_cbar_extend']

        if i%2 == 0:
            cbar = m.colorbar(cs, location='left',extend = extend )
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.yaxis.set_ticks_position('left')
        else:
            cbar = m.colorbar(cs,location='right',extend=extend)
        cbar.set_label(units)

        plt.title(Var_Dict[pannel+'_label'], fontsize=9)
            
    fig.suptitle(Var_Dict['plot_title'], fontsize=14, fontweight='bold')
    fig.text(0.5, 0.92, subtitle, ha='center',size='medium')
    figout = os.path.join(Main['outpath'],outname)

    plt.savefig(figout, format=Main['outformat'], dpi=Main['dpi'])
    print 'done with',figout
    return

if __name__ == "__main__":
    main()    
