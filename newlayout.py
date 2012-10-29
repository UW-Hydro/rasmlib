#!/usr/bin/env python

import argparse, time
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import pdb

import pickle
vari=pickle.load(open('vab.p','rb'))
ncdata=pickle.load(open('ncdata.p','rb'))


def main(ncdata,vari):

        vabs = vari['vabs']
        units = vari['units']
        pathout = vari['pathout']
        dayin = vari['dayin']
        colormp = vari['colormp']
        subtitles = vari['subtitles']
        title = vari['title']
        outformat = vari['outformat']

        print('the output directory is: '+pathout)

        create_plot(ncdata,vabs,units,pathout,dayin,colormp,subtitles,title,outformat)



def create_plot(ncdata,vabs,units,pathout,dayin,colormp,subtitles,title,outformat):

        i=0
        for varb in vabs:
            data = ncdata[varb][:]
            unit = units[i]
            i = i+1
            lon = ncdata['lon'][:]
            lat = ncdata['lat'][:]
#            nlat = np.ma.compressed(lat)
#            latmi = min(nlat)

            plt.subplot(2,2,i)
            m = Basemap(projection='npstere',boundinglat=20,lon_0=-114,\
                        lat_ts=80.5,resolution='i')
            m.drawcoastlines()
            m.drawparallels(np.arange(20,91,20))
            m.drawmeridians(np.arange(-180,181,20))



            xi,yi = m(lon,lat)
            cm = matplotlib.cm.get_cmap(colormp)
            cs = m.pcolor(xi,yi,data,cmap=cm)
            plt.title(subtitles[i-1])
            cbar = m.colorbar(cs,location='right')
            cbar.set_label(unit) 

        plt.suptitle(title)
        outnamme = title+'_'+dayin+'.'+outformat
        figout=os.path.join(pathout,outnamme)
        plt.savefig(figout,format=outformat)
        return

if __name__=="__main__":
        main(ncdata,vari)

