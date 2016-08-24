'''
This module studies radiation at a global scale.
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import calendar as cal
import time
from ent_datasets import months_in_year

# Load datasets
grid = np.load('datasets/earth_grid.npy')
lat, lon = grid[:,0,0], grid[0,:,1]

wvlen = np.load('datasets/wvlen.npy')           # nm
wvlen_lres = np.load('datasets/wvlen_lres.npy') # nm
wvnum = np.load('datasets/wvnum.npy')           # cm^-1
wvlen_num = 1.0e7/wvnum                         # nm

titles = ['SW (low res)', 'LW', 'Clear sky LW']

def loadflux(month, re='r'):
    '''
    Loads SW and LW radiation datasets which correspond to given month.
    
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    Parameter re takes the values 'r' and 'e' for radiation and entropy 
    fluxes respectively.
    
    SW radiation is corrected by a factor of 1000 if not lres.
    SW radiation units: uW cm^-2 sr^-1 nm^-1
    LW radiation units: W cm^-2 sr^-1 cm
    '''
    swl = np.load('datasets/flux/lres_sw%s%s.npy' % (re,month))
    lw = np.load('datasets/flux/lw%s%s.npy' % (re,month))
    lwc = np.load('datasets/flux/clr_lw%s%s.npy' % (re,month))
    return swl, lw, lwc

# Plot functions
def plot_flux(month, flux, re='r'):
    '''
    This function plots global maps of given SW and LW fluxes.
    
    Parameter flux is an array of fluxes.
    '''
    s, u = is_re(re)
    m = calendar(month)
    xx, yy, flux = shift_grid(flux)
    
    fig = plt.figure()
    fig.suptitle("Global "+s+" flux for "+m[1]+" "+m[0], fontsize=16)
    for i in range(len(flux)):
        ax = fig.add_subplot(2,2,i+1)
        ax.set_title(titles[i]+" flux", fontsize=14)
        m = Basemap()
        m.drawmapboundary()
        m.drawcoastlines()
        m.contourf(xx,yy,flux[i])
        m.colorbar(location='bottom', label='Flux '+u)
    plt.tight_layout()
    plt.subplots_adjust(top=0.925)

def plot_diff(month1, month2, flux1, re='r', info=False):
    '''
    Parameter month2 corresponds to the month, month1 (with flux1) is
    compared with.
    
    This function plots the difference of fluxes between two months.
    If info=True, the difference is returned and not plotted.
    '''
    s, u = is_re(re)
    m1, m2 = calendar(month1), calendar(month2)
    flux2 = np.array(loadflux(month2, re))
    
    fig = plt.figure()
    fig.suptitle("Difference of "+s+" flux between "+m1[1]+" "
                 +m1[0]+" and "+m2[1]+" "+m2[0], fontsize=16)
    
    xx, yy, flux1 = shift_grid(flux1)
    xx, yy, flux2 = shift_grid(flux2)
    fdiff = flux1-flux2
    if info: return xx, yy, fdiff
    
    for i in range(len(fdiff)):
        ax = fig.add_subplot(2,2,i+1)
        ax.set_title(titles[i]+" flux difference", fontsize=14)
        m = Basemap()
        m.drawmapboundary()
        m.drawcoastlines()
        m.contourf(xx,yy,fdiff[i])
        m.colorbar(location='bottom', label='Flux '+u)
    plt.tight_layout()
    plt.subplots_adjust(top=0.925)

# Utility functions
def find_nans(arr):
    nans = np.where(np.isnan(arr))
    nans = map(tuple,np.array(nans).T)
    return nans

def sumup(arr):
    return {'shape':arr.shape, 'max':np.nanmax(arr), 'min':np.nanmin(arr),
            'nans':np.where(np.isnan(arr))[0].size}

def calendar(month):
    y, m = month[:2], int(month[2:])
    if y=='00': y = '2000'
    elif y=='99': y = '2099'
    m = cal.month_name[m]
    return y, m

def is_re(re):
    if re=='r': s, u = 'radiation', '($W\ m^{-2}\ sr^{-1}$)'
    elif re=='e': s, u = 'entropy', '($mW\ m^{-2}\ sr^{-1}\ K^{-1}$)'
    return s, u

def shift_grid(flux):
    for i in range(len(flux)):
        flux[i], lons = shiftgrid(180.0, flux[i], lon, start=False)
    xx, yy = np.meshgrid(lons, lat)
    return xx, yy, flux

# Main functions
def analyse_month(month, re='r', diff=None, gmap=True, info=False):
    '''
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01 = January).
    
    This function plots SW & LW radiation and entropy fluxes across the 
    globe. It can also provide useful information about the datasets
    (e.g. extreme values)
    '''
    flux = np.array(loadflux(month, re))
    if gmap:
        plot_flux(month, flux, re)
        if diff is not None:
            for m in diff:
                plot_diff(month, m, flux, re, info=False)
    
    if info:
        return sumup(sw_eflux), sumup(lw_eflux)

analyse_month('0009', 'e', diff=['0008','0010','9909'])
plt.show()


