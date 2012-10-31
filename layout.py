#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pylab import *
import pdb
import os

import pickle
a=pickle.load(open('ncdata.p','rb'))
b=pickle.load(open('vari.p','rb'))

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
        
        locs = ['left','right','left','right']
        fig = plt.figure()
        b = 0.04;
        c = 0;
        b = 0.5*b;
        c = 0.5*c;
        xxs = [0.15+b, 0.5-b, 0.15+b, 0.5-b]
        yys = [0.525-c, 0.525-c, 0.125+c, 0.125+c]
        subwidth = 0.35
        subheight = 0.33
        i=0
        while i<=3:
            varb = vabs[i]
            data = ncdata[varb][:]
            unit = units[i]
            axes([xxs[i], yys[i], subwidth, subheight])
            i = i+1
            lon = ncdata['lon'][:]
            lat = ncdata['lat'][:]

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
                    cbar = m.colorbar(cs,location='left')
                    cbar.ax.yaxis.set_label_position('left')
                    cbar.ax.yaxis.set_ticks_position('left')
            else:
                    cbar = m.colorbar(cs,location='right')
            cbar.set_label(unit)

        fig.suptitle(title,fontsize=14, fontweight='bold')
        fig.text(0.5, 0.90, subtitle, ha='center',size='medium')
        outname = outnamme+'.'+outformat
        figout = os.path.join(pathout,outname)
        plt.show()
        plt.savefig(figout,format=outformat)
        return
if __name__ == "__main__":
        main(a,b)
