#!/usr/bin/env python

import argparse
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("inputp", help="the directory that stores the files")
        parser.add_argument("date", help="the date that you want such as 1999-06-26")
        parser.add_argument("outputp", help="the directory that you want to save the output picture")
        args = parser.parse_args()
        
        pathin=args.inputp
        dayin=args.date
        pathout=args.outputp
        filepath=pathin+'r32RB1a.vic.h.'+dayin+'-00000.nc'
        print(filepath)
        print(pathout)
        if error_check(pathout,filepath):
                return
        nc=read_netcdf(filepath)
        create_plot(nc,pathout,dayin)


def error_check(pathout,filepath):

        if not(os.path.isdir(pathout)):
                print('the output directory do not exist')
                return(1)
        if not(os.path.exists(filepath)):
                print('the netcdf file do not exist in the directory you typed')
                return(1)
        return

def read_netcdf(filepath):

        nc=Dataset(filepath)
        return(nc)

    
def create_plot(nc,pathout,dayin):
        i=0
        vabs=['Precipitation','Swq','Swnet','Tair']
        for varb in vabs:
            rawdata=nc.variables[varb]
            data=rawdata[0][:]
            ndata=np.ma.compressed(data)
            unit=str(rawdata.units)
            i=i+1
            lon=nc.variables['longitude'][:]
            lat=nc.variables['latitude'][:]
            nlat=np.ma.compressed(lat)
            latmi=min(nlat)

            plt.subplot(2,2,i)
            m=Basemap(projection='nplaea',boundinglat=latmi,lon_0=160,
                      resolution='l')
            m.drawcoastlines()
            m.drawparallels(np.arange(20,91,20))
            m.drawmeridians(np.arange(-180,181,20))



            xi,yi=m(lon,lat)
            cs=m.contourf(xi,yi,data)
            plt.title(varb)
            cbar=m.colorbar(cs,location='right')
            cbar.set_label(unit) # xunhuan danwei - don't use chinese comments

            
            
        tit='RASM Data_'
        titwithdate=tit+dayin
        plt.suptitle(titwithdate)
        figout=pathout+dayin+'.png'
        plt.savefig(figout,format='png')
        return

if __name__=="__main__":
	main()

