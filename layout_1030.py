#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import os


def main(ncdata,vari):

        vabs = vari['vabs']
        units = vari['units']
        projection = vari['projection']
        lorange = vari['lorange']
        larange = vari['larange']
        colormp = vari['colormp']
        labels = vari['labels']
        subtitle = vari['subtitle']
        title = vari['title']
        pathout = vari['pathout']
        outnamme = vari['outnamme']
        outformat = vari['outformat']

        print('the output directory is: '+pathout)

        create_plot(ncdata,vabs,units,projection,lorange,larange,colormp,\
                    labels,subtitle,title,pathout,outnamme,outformat)



def create_plot(ncdata,vabs,units,projection,lorange,larange,colormp,labels,\
                subtitle,title,pathout,outnamme,outformat):
        
#        rc('text', usetex=True)
#        rc('font', family='serif')
        locs = ['left','right','left','right']
        fig = plt.figure()
        plt.subplots_adjust(top=0.85, right=0.87)

        i=0
        while i<=3:
            varb = vabs[i]
            data = ncdata[varb][:]
            unit = units[i]
            i = i+1
            lon = ncdata['lon'][:]
            lat = ncdata['lat'][:]
#            nlat = np.ma.compressed(lat)
#            latmi = min(nlat)

            plt.subplot(2,2,i)
            m = Basemap(projection=projection,boundinglat=20,lon_0=-114,\
                        lat_ts=80.5,resolution='i')
            m.drawcoastlines()
            m.drawparallels(np.arange(*lorange))
            m.drawmeridians(np.arange(*larange))

            xi,yi = m(lon,lat)
            cm = matplotlib.cm.get_cmap(colormp)
            cs = m.pcolor(xi,yi,data,cmap=cm)
            plt.title(labels[i-1])
            if i%2 == 1:
                    cbar = m.colorbar(cs,location='left',pad=0.58)
            else:
                    cbar = m.colorbar(cs,location='right',pad=0.25)
            cbar.set_label(unit) 

#        suptitle = '\\textbf{'+title+'}'
#        subtitle = '{\\small \\textrm{'+subtitle+'}}'
        fig.suptitle(title,fontsize=14, fontweight='bold')
        fig.text(0.5, 0.90, subtitle, ha='center',size='medium')
        outname = outnamme+'.'+outformat
        figout = os.path.join(pathout,outname)
        plt.savefig(figout,format=outformat)
        return

