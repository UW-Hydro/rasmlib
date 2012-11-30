import numpy as np
import matplotlib
from matplotlib.cm import *

def cmap_discretize(cmap, N): 

     cmap = get_cmap(cmap)
     colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
     colors_rgba = cmap(colors_i)
     indices = np.linspace(0, 1., N+1)
     cdict = {}
     for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])\
                     for i in xrange(N+1) ]

     return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


