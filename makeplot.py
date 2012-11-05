#!/usr/bin/env python

from netCDF4 import Dataset
from layout  import create_plot
import argparse
import numpy.ma as ma

def main():
        variname = 'Air Temperature'
        varisavename = 'air_temperature'
        rvari = 'Tair'
        evari = 't2m'

        args = processcml()
        ncr1 = Dataset(args.r1)
        ncr2 = Dataset(args.r2)
        ncr3 = Dataset(args.r3)
        ncera = Dataset(args.era)
        ncdatam = {}
        ncdatam['lon'] = ncr1.variables['longitude'][:]
        ncdatam['lat'] = ncr1.variables['latitude'][:]
        ncdatam['r1'] = ncr1.variables[rvari][0][:]
        ncdatam['r2'] = ncr2.variables[rvari][0][:]
        ncdatam['r3'] = ncr3.variables[rvari][0][:]
        ncdatam['era'] = ncera.variables[evari][0][:]-272.15 # K to C units
        ncdatad = {}
        ncdatad['lon'] = ncdatam['lon']
        ncdatad['lat'] = ncdatam['lat']
        ncdatad['r1'] = ncdatam['era']-ncdatam['r1']
        ncdatad['r2'] = ncdatam['era']-ncdatam['r2']
        ncdatad['r3'] = ncdatam['era']-ncdatam['r3']
        ncdatad['era'] = ncdatam['era']

        unit = ncr1.variables[rvari].units


        vabs = ['r1','r2','r3','era']
        units = [unit for i in range(4)]
        projection = 'npstere'
        lorange = [-80, 81, 20]
        larange = [-180, 181, 20]
        colormp = 'jet'
        labels = ['R30RB1g','R33RBVIC60','R33RBVIC70','ERA']
        subtitle = '1991-01'
        titlem = variname+' mean'
        titled = variname+' anomaly'
        pathout = args.figpath
        outnammem = varisavename+'_'+subtitle+'_mean'
        outnammed = varisavename+'_'+subtitle+'_anomaly'
        outformat = 'png'
        projection_parameters = {'projection': 'npstere', 'boundinglat': 35,\
                                 'lon_0': -114, 'lat_ts': 80.5}

        maxm = max([ma.MaskedArray.max(ncdatam[vab]) for vab in vabs])
        minm = min([ma.MaskedArray.min(ncdatam[vab]) for vab in vabs])
        maxd = max([ma.MaskedArray.max(ncdatad[vab]) for vab in vabs])
        mind = min([ma.MaskedArray.min(ncdatad[vab]) for vab in vabs])

        print('the output directory is: '+pathout)

        create_plot(ncdatam,vabs,units,lorange,larange,colormp,minm,maxm,\
                    labels,subtitle,titlem,pathout,outnammem,outformat,\
                    projection_parameters)

        create_plot(ncdatad,vabs,units,lorange,larange,colormp,mind,maxd,\
                    labels,subtitle,titled,pathout,outnammed,outformat,\
                    projection_parameters)

        return


def processcml():

    parser = argparse.ArgumentParser(description='Plot RASM and ERA interim albedo')
    parser.add_argument('--r1', required=True, metavar='<RASM VIC netcdf file>',
                        help='RASM VIC netCDF file plotted in upper left')
    parser.add_argument('--r2', required=True, metavar='<RASM VIC netcdf file>',
                        help='RASM VIC netCDF file plotted in upper right')
    parser.add_argument('--r3', required=True, metavar='<RASM VIC netcdf file>',
                        help='RASM VIC netCDF file plotted in lower left')
    parser.add_argument('--era', required=True, metavar='<ERA interim netcdf file>',
                        help='ERA Interim file for same period as RASM VIC files')
    parser.add_argument('--figpath', required=True, metavar='<path for output files>',
                        help='Directory for output files')
    return parser.parse_args()

