'''
This module studies radiation at a global scale
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap

# Load datasets
# SW rad is corrected by a factor of 1000
grid = np.load('datasets/earth_grid.npy')

wvlen = np.load('datasets/wvlen.npy') # nm
wvnum = np.load('datasets/wvnum.npy') # cm^-1
wvlen_num = 1.0e7/wvnum               # nm

sw0001 = 1000*np.load('datasets/sw0001.npy') # uW cm^-2 sr^-1 nm^-1
lw0001 = np.load('datasets/lw0001.npy')      # W cm^-2 sr^-1 cm

# Physical constants
c=2.998e8   # m s^-1
kb=1.38e-23 # kg m^2 s^-2 K^-1
h=6.626e-34 # J s = kg m^2 s^-1

# Conversion functions
def radtoent(rad,radtype='sw'):
    '''
    Converts array of radiation intensity to array of entropy.
    
    Needs radiation in W m^-2 sr^-1 m.
    Returns entropy in mW m^-2 sr^-1 K^-1 nm^-1.
    '''
    if radtype=='sw':
        wvl = wvlen
        iconst = 1.0e-11*wvl*wvl
        iconst = iconst[:,None,None]
    elif radtype=='lw':
        wvl = wvlen_num
        iconst = 1.0e2
    wvn = 1.0e9/wvl
    econst = 1.0e-6*wvn*wvn
    econst, wvn = econst[:,None,None], wvn[:,None,None]
    
    intens = iconst*rad
    y = intens/(2*h*c*c*wvn**3)
    ent = np.where((y!=0.0).any(),
          econst*2*kb*c*wvn*wvn*((1+y)*np.log(1+y)-y*np.log(y)),0.0)
    
    return ent

def radtorad(rad,radtype='sw'):
    '''
    Converts radiation to units W m^-2 sr^-1 nm^-1
    '''
    if radtype=='sw':
        wvln = wvlen
        rconst = 1.0e-2
    elif radtype=='lw':
        wvln = wvnum
        rconst = 1.0e-3*wvnum*wvnum
        rconst = rconst[:,None,None]
    
    return rconst*rad

m = Basemap()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcoastlines()
plt.show()








