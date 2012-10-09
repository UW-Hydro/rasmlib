#!/usr/bin/env python

import argparse, time
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("inputp", help="the directory that stores the files")
        parser.add_argument("case", help="the case that you want such as r32RB1a")
        parser.add_argument("date", help="the date that you want such as 1999-06-26")
        parser.add_argument("outputp", help="the directory that you want to save the output picture")
        args = parser.parse_args()
        
        pathin=args.inputp
        dayin=args.date
        pathout=args.outputp
        case=args.case
        namme=case+'.vic.h.'+dayin+'-00000.nc'
        filepath=os.path.join(pathin,namme)
        print(filepath)
        print(pathout)
        if error_check(pathout,filepath,dayin):
                return
#        if datecheck(dayin):
#                print('the date you typed is not valid')
#                return
        nc=read_netcdf(filepath)
        create_plot(nc,pathout,dayin)


def error_check(pathout,filepath,dayin):

        try:
            valid_date=time.strptime(dayin,'%Y-%m-%d')
        except ValueError:
            print('the date you gave is wrong')
            return(1)
        try:
            nc=Dataset(filepath)
        except RuntimeError:
	    print('the file does not exist in the directory you typed')
            return(1)
        if not(os.path.isdir(pathout)):
            print('the output directory does not exist')
            return(1)
        return(0)

def datecheck(date):
    year=int(date[0:4])
    mon=int(date[5:7])
    dat=int(date[8:10])
    if mon<1 or mon > 12:
       return(1)
    if dat>31 or dat < 1:
       return(1)
    if (mon-4)*(mon-6)*(mon-9)*(mon-11)==0:
       if dat>30:
          return(1)
    if mon==2:
       if dat>29:
          return(1)
       if dat==29:
          if year%4 != 0:
             return(1)
          if year%100 == 0 and year%400 != 0:
             return(1)
    return(0)

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
            m=Basemap(projection='nplaea',boundinglat=latmi,lon_0=160,\
                      resolution='l')
            m.drawcoastlines()
            m.drawparallels(np.arange(20,91,20))
            m.drawmeridians(np.arange(-180,181,20))



            xi,yi=m(lon,lat)
            cs=m.contourf(xi,yi,data)
            plt.title(varb)
            cbar=m.colorbar(cs,location='right')
            cbar.set_label(unit) 

            
            
        tit='RASM Data_'
        titwithdate=tit+dayin
        plt.suptitle(titwithdate)
        outnamme=dayin+'.png'
        figout=os.path.join(pathout,outnamme)
        plt.savefig(figout,format='png')
        return

if __name__=="__main__":
	main()

