"""
rasmlib plot utilities
"""

def sub_plot_pcolor(lons, lats, data, title=None, cmap=cm.jet,
                    vmin=None, vmax=None, cbar=True, cbar_location='bottom',
                    units=None):
    

    if vmin is None:
        vmin=data.min()
    if vmax is None:
        vmax=data.max()

    projection = {'urcrnrlat': 27.511827255753555, 
                  'urcrnrlon': 16.90845094934209, 
                  'llcrnrlat': 16.534986367884521, 
                  'llcrnrlon': 189.2229322311162, 
                  'projection': 'lcc', 
                  'rsphere': 6371200.0, 
                  'lon_0': -114, 
                  'lat_0': 90}
        
    m = Basemap(**projection)
    xi,yi = m(np.squeeze(lons),np.squeeze(lats))
    sp = m.pcolormesh(xi, yi, np.squeeze(data), 
                      vmin=vmin, vmax=vmax, cmap=cmap)
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawcoastlines(color='k',linewidth=0.25);
    
    if title:
        plt.title(title,size=13)
    if cbar:
        cbar = m.colorbar(location=cbar_location)
    cbar.set_label(units)
    
    return sp

def cmap_discretize(cmap, N=10):

     cmap = cm.get_cmap(eval(cmap))
     colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
     colors_rgba = cmap(colors_i)
     indices = np.linspace(0, 1., N+1)
     cdict = {}
     for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1)]

     return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
    
def plot4(lons, lats, pannels, variable, titles=None, 
          vmin=None, vmax=None, units=None, cmap=None, 
          amin=None, amax=None, amap=None,
          outdir=None, y0=None, y1=None, mask=None, relative=True):
    """
    Make 4 pannel plots for absolute values and anomaly fields
    lons - longitudes (2d array)
    lats - latitudes (2d array)
    pannels - dictionary of data pannels[pannel_number][variable][season, y, x]
    titles - optional names for each pannel
    vmin - 
    vmax - 
    """
    
    if titles is None:
        titles = pannels.keys()
    
    for i, season in enumerate(seasons):
        if y0 is not None and y1 is not None:
            title = '%s-%s %s %s' %(y0, y1, season, variable)
        else:
            title = '%s-%s' %(season, variable)
        fig = plt.figure(figsize=(10,9))
        for j in xrange(1, 5):
            ax = fig.add_subplot(2,2,j)
            sp = sub_plot_pcolor(lons, lats, 
                                 np.ma.masked_where(mask==0, pannels[j][variable].sel(season=season)), 
                                 title=titles[j], 
                                 vmin=vmin, vmax=vmax, units=units , cmap=cmap)
        fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
        fig.tight_layout()
        
        if outdir:
            plt.savefig(os.path.join(outdir, '%s-%s-%s.png' %(variable, season, run)), 
                        dpi=300, format='png', bbox_inches='tight')

        fig = plt.figure(figsize=(10,9))
        for j in xrange(1, 5):
            ax = fig.add_subplot(2,2,j)
            if j == 4:
                sp = sub_plot_pcolor(lons, lats,
                                     np.ma.masked_where(mask==0, pannels[j][variable].sel(season=season)), 
                                     title=titles[j], vmin=vmin, vmax=vmax, 
                                     units=units , cmap=cmap)
            else:
                diff = pannels[j][variable].sel(season=season).values \
                    - pannels[4][variable].sel(season=season).values
                sp = sub_plot_pcolor(lons, lats, 
                                     np.ma.masked_where(mask==0, diff), 
                                     title='%s-%s' %(titles[j], titles[4]), 
                                     vmin=amin, vmax=amax, units=units,
                                     cmap=amap)
        fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
        fig.tight_layout()
        
        if outdir:
            plt.savefig(os.path.join(outdir, '%s-anom-%s-%s.png' % (variable, season, run)), 
                        dpi=300, format='png', bbox_inches='tight')

        if relative:
            fig = plt.figure(figsize=(10,9))
            for j in xrange(1, 5):
                ax = fig.add_subplot(2,2,j)
                if j == 4:
                    sp = sub_plot_pcolor(lons, lats, np.ma.masked_where(mask==0,
                                         pannels[j][variable].sel(season=season)), 
                                         title=titles[j], vmin=vmin, vmax=vmax, 
                                         units=units , cmap=cmap)
                else:
                    diff = 100*(pannels[j][variable].sel(season=season).values- pannels[4][variable].sel(season=season).values)/pannels[4][variable].sel(season=season).values
                    sp = sub_plot_pcolor(lons, lats, np.ma.masked_where(mask==0, diff), 
                                         title='(%s-%s)/%s' %(titles[j], titles[4], titles[4]), 
                                         vmin=-100, vmax=100, units='%' , cmap=amap)
            fig.suptitle(title, fontsize=14, fontweight='bold', y=1.03)
            fig.tight_layout()

            if outdir:
                plt.savefig(os.path.join(outdir, '%s-relative-anom-%s-%s.png' %(variable, season, run)), 
                            dpi=300, format='png', bbox_inches='tight')
