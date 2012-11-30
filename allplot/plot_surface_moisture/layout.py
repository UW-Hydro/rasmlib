#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from pylab import *
import pdb
import os
from matplotlib.cm import *
from cut_colmap import cmap_discretize

def create_plot(kind,ncdata,vabs,units,lorange,larange,colormp,minv,maxv,labels,\
                subtitle,title,pathout,outnamme,outformat,projection_parameters):
        
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

            m = Basemap(**projection_parameters)
            m.drawcoastlines(linewidth=0.25)
            m.drawparallels(np.arange(*lorange))
            m.drawmeridians(np.arange(*larange))

            xi,yi = m(lon,lat)

            if kind == 1 or i ==4:
                 cm = get_cmap(colormp)
                 cs = m.pcolor(xi, yi, data, cmap=cm, vmin=minv, vmax=maxv)
            else:
                 cm = get_cmap('RdYlBu_r')
                 vv = np.linspace(-2, 10, -10, endpoint=True)
                 cm = cmap_discretize(cm, 10)
                 wd = (maxv-minv)*0.25
                 cs = m.pcolor(xi, yi, data, cmap=cm, vmin=-wd, vmax=wd)

            plt.title(labels[i-1], fontsize=9)

            if i%2 == 1:
                    cbar = m.colorbar(cs, location='left')
                    cbar.ax.yaxis.set_label_position('left')
                    cbar.ax.yaxis.set_ticks_position('left')
            else:
                    cbar = m.colorbar(cs,location='right')
            cbar.set_label(unit, fontsize=9)
            for t in cbar.ax.get_yticklabels():
                    t.set_fontsize(9)

        fig.suptitle(title, fontsize=14, fontweight='bold')
        fig.text(0.5, 0.90, subtitle, ha='center', fontsize=11)# size='medium')
        outname = outnamme+'.'+outformat
        figout = os.path.join(pathout, outname)
        plt.savefig(figout, format=outformat, dpi=200)
        return


