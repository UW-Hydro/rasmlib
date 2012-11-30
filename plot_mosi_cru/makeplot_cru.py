#!/usr/bin/env python

from netCDF4 import Dataset
from layout  import create_plot
import argparse
import numpy.ma as ma
import pdb
import os

def main(variname, varisavename, rvari, evari):

        args = processcml()
        ncr1 = Dataset(args.r1)
        ncr2 = Dataset(args.r2)
        ncr3 = Dataset(args.r3)
        ncera = Dataset(args.era)
        timee = args.time
        
        if os.path.isfile(args.r1) and os.path.isfile(args.r2) and\
           os.path.isfile(args.r3) and os.path.isfile(args.era):
           print('\n ok')
        else: 
           print('\n input directory is wrong')
           return
        if not(os.path.exists(args.figpath)):
           print('\n output directory is wrong')

        if evari == 't2m':
           cvari = 'tmp'
        if evari == 'tp':
           cvari = 'prcp';

        ncdatam = {}
        ncdatam['lon'] = ncr1.variables['longitude'][:]
        ncdatam['lat'] = ncr1.variables['latitude'][:]
        ncdatam['r1'] = ncr1.variables[rvari][0][:]
        ncdatam['r2'] = ncr2.variables[rvari][0][:]
        ncdatam['r3'] = ncr3.variables[evari][:]  # this now actually era
        ncdatam['era'] = ncera.variables[cvari][:] # this now actually cru

        if evari == 't2m':
           ncdatam['r3'] = ncdatam['r3']-273.15
        if evari == 'tp':
           ncdatam['r3'] = ncdatam['r3']*1000

        ncdatad = {}
        ncdatad['lon'] = ncdatam['lon']
        ncdatad['lat'] = ncdatam['lat']
        ncdatad['r1'] = -ncdatam['era']+ncdatam['r1']
        ncdatad['r2'] = -ncdatam['era']+ncdatam['r2']
        ncdatad['r3'] = -ncdatam['era']+ncdatam['r3']
        ncdatad['era'] = ncdatam['era']

        unit = ncr1.variables[rvari].units


        vabs = ['r1','r2','r3','era']
        units = [unit for i in range(4)]
        projection = 'npstere'
        lorange = [-80, 81, 20]
        larange = [-180, 181, 20]
        colormp = 'jet'
        labels = ['VIC60','VIC70','ERA','CRU']
        subtitle = timee
        titlem = variname+' mean'
        titled = variname+' anomaly'
        pathout = args.figpath
        outnammem = varisavename+'_'+subtitle+'_mean'
        outnammed = varisavename+'_'+subtitle+'_anomaly'
        outformat = 'png'
        projection_parameters = {'projection': 'npstere', 'boundinglat': 53,\
                                 'lon_0': -114, 'lat_ts': 80.5}

        maxm = max([ ma.MaskedArray.max(ncdatam[vab]) for vab in vabs[:3] ])
        minm = min([ ma.MaskedArray.min(ncdatam[vab]) for vab in vabs[:3] ])
        mxmi = maxm - minm
        maxd = 0.25*mxmi
        mind = -0.25*mxmi

        if rvari == 'Swq'  and evari == 'sd':
           maxm = 300
        if evari == 'tp':
           maxm = 5
           

        print('the output directory is: '+pathout)

        create_plot(1,ncdatam,vabs,units,lorange,larange,colormp,minm,maxm,\
                    labels,subtitle,titlem,pathout,outnammem,outformat,\
                    projection_parameters)

        create_plot(2,ncdatad,vabs,units,lorange,larange,colormp,minm,maxm,\
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
    parser.add_argument('--time', required=True, metavar='<time consist with type>',
                        help='the time period of the plot')    

    return parser.parse_args()

