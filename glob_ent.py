'''
This module studies radiation at a global scale.
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import calendar as cal
from ent_datasets import months_in_year

import time

# Load datasets
grid = np.load('datasets/earth_grid.npy')
lat, lon = grid[:,0,0], grid[0,:,1]

wvlen = np.load('datasets/wvlen.npy')           # nm
wvlen_lres = np.load('datasets/wvlen_lres.npy') # nm
wvnum = np.load('datasets/wvnum.npy')           # cm^-1
wvlen_num = 1.0e7/wvnum                         # nm

def loadflux(month, re='r'):
    '''
    Loads SW and LW radiation datasets which correspond to given month.
    
    SW radiation is corrected by a factor of 1000 if not lres.
    SW radiation units: uW cm^-2 sr^-1 nm^-1
    LW radiation units: W cm^-2 sr^-1 cm
    
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    Parameter re takes the values 'r' and 'e' for radiation and entropy 
    fluxes respectively.
    '''
    swl = np.load('datasets/flux/lres_sw%s%s.npy' % (re,month))
    lw = np.load('datasets/flux/lw%s%s.npy' % (re,month))
    lwc = np.load('datasets/flux/clr_lw%s%s.npy' % (re,month))
    return swl, lw, lwc

titles = ['SW (low res)', 'LW', 'Clear sky LW']
rad = {'swl':0, 'lw':1, 'lwc':2}

# Plot functions
def plot_map(xx, yy, data, u):
    '''
    Plots a global map of the data provided.
    
    Parameters correspond to londitude and latitude, data to be plotted
    and units.
    '''
    m = Basemap()
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,data)
    m.colorbar(location='bottom', label='Flux '+u)

def plot_flux(month, flux, re='r'):
    '''
    Plots global maps of given SW and LW fluxes.
    
    Parameter flux is an array of fluxes.
    '''
    s, u = is_re(re)
    m = calendar(month)
    xx, yy, shflux = shift_grid(flux)
    
    fig = plt.figure()
    fig.suptitle("Global "+s+" flux for "+m[1]+" "+m[0], fontsize=16)
    for i in range(len(shflux)):
        ax = fig.add_subplot(2,2,i+2)
        ax.set_title(titles[i]+" flux", fontsize=14)
        plot_map(xx, yy, shflux[i], u)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.925)

def plot_diff(month1, month2, flux1, re='r', info=False):
    '''
    Plots the difference of fluxes between two months.
    If info=True, the difference is returned and not plotted.
    
    Parameter month2 corresponds to the month, month1 (with flux1) is
    compared with.
    '''
    s, u = is_re(re)
    m1, m2 = calendar(month1), calendar(month2)
    flux2 = np.array(loadflux(month2, re))
    
    fig = plt.figure()
    fig.suptitle("Difference of "+s+" flux between "+m1[1]+" "
                 +m1[0]+" and "+m2[1]+" "+m2[0], fontsize=16)
    
    fdiff = flux1-flux2
    xx, yy, shfdiff = shift_grid(fdiff)
    if info: return xx, yy, shfdiff
    
    for i in range(len(fdiff)):
        ax = fig.add_subplot(2,2,i+2)
        ax.set_title(titles[i]+" flux difference", fontsize=14)
        plot_map(xx, yy, shfdiff[i], u)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.925)

# Utility functions
def find_nans(arr):
    '''
    Returns the indices of all nans in an array arr.
    '''
    nans = np.where(np.isnan(arr))
    nans = map(tuple,np.array(nans).T)
    return nans

def sumup(arr):
    '''
    Returns useful information about array arr (eg. extreme values)s.
    '''
    return {'shape':arr.shape, 'max':np.nanmax(arr), 'min':np.nanmin(arr),
            'nans':np.where(np.isnan(arr))[0].size}

def calendar(month):
    '''
    Converts a string of the form yymm (y for year, m for month) to yyyy
    and the name of the corresponding month m (eg. input '0001' returns
    output y, m = '2000', 'January'.
    '''
    y, m = month[:2], int(month[2:])
    if y=='00': y = '2000'
    elif y=='99': y = '2099'
    m = cal.month_name[m]
    return y, m

def neighbour_months(month):
    '''
    Returns the previous and next month as well as the same month in the other
    year (2000 or 2099).
    '''
    if month[:2]=='00': yr, notyr = 2000, 2099
    elif month[:2]=='99': yr, notyr = 2099, 2000
    months = months_in_year(yr)+months_in_year(notyr)
    ind = months.index(month)
    if ind==0: inds = [12,1]
    elif ind==11: inds = [23,10]
    else: inds = [ind+12,ind-1,ind+1]
    neighbours = []
    for ind in inds:
        neighbours.append(months[ind])
    return neighbours

def is_re(re):
    '''
    Returns 'radiation' and 'entropy' to be used in titles and plot labels
    for re = 'r' and 'e' respectively, with the corresponding flux units.
    '''
    if re=='r': s, u = 'radiation', '($W\ m^{-2}\ sr^{-1}$)'
    elif re=='e': s, u = 'entropy', '($mW\ m^{-2}\ sr^{-1}\ K^{-1}$)'
    return s, u

def shift_grid(flux):
    '''
    Shifts longitude and data arrays to match a map centred at
    (lat, lon) = (0, 0).
    '''
    shift_flux = np.empty(flux.shape)
    for i in range(len(flux)):
        shift_flux[i], lons = shiftgrid(180.0, flux[i], lon, start=False)
    xx, yy = np.meshgrid(lons, lat)
    return xx, yy, shift_flux

# Main functions
def analyse_month(month, re='r', diff=None, gmap=True, info=False):
    '''
    Plots SW & LW radiation and entropy fluxes across the 
    globe. It can also call function sumup for the flux data of the month.
    
    Parameter diff is a list of months to be compared with month, gmap and info
    plot the relevant maps and provide sumup information if True respectively.
    '''
    flux = np.array(loadflux(month, re))
    print flux.shape
    if gmap:
        plot_flux(month, flux, re)
        if diff is not None:
            for m in diff:
                plot_diff(month, m, flux, re, info=False)
    
    if info:
        temp = []
        for f in flux: temp.append(sumup(f))
        return temp

def spec_analysis(month, re='r', radtype='swl', gmap=True, info=False):
    '''
    Plots monthly comparisons of specific type of flux across the globe
    (SW, LW, Clear sky LW).
    
    Parameter radtype ('swl', 'lw', 'lwc') determines the type of flux.
    '''
    flux = np.array([loadflux(month, re)[rad[radtype]]])
    s, u = is_re(re)
    neighbours = neighbour_months(month)
    m = calendar(month)
    
    fig = plt.figure()
    fig.suptitle("Comparison of "+titles[rad[radtype]]+" "+s+" flux for "
    +m[1]+" "+m[0], fontsize=16)
    
    ax = fig.add_subplot(2,2,1)
    ax.set_title("Flux in "+m[1]+" "+m[0], fontsize=14)
    xx, yy, shflux = shift_grid(flux)
    if gmap: plot_map(xx, yy, shflux[0], u)
    
    fdiffs = []
    for n in range(len(neighbours)):
        ax = fig.add_subplot(2,2,n+2)
        m2 = calendar(neighbours[n])
        ax.set_title("Difference with "+m2[1]+" "+m2[0], fontsize=14)
        flux2 = np.array([loadflux(neighbours[n], re)[rad[radtype]]])
        fdiff = flux - flux2
        xx, yy, fdiff = shift_grid(fdiff)
        fdiffs.append((xx, yy, fdiff[0]))
        if gmap: plot_map(xx, yy, fdiff[0], u)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.925)
    if info: return fdiffs

analyse_month('9902', re='r', diff=None)
#plot_diff('9907','9908', loadflux('9907', 'e'), re='e')
#spec_analysis('9902', 'r', 'lw')
plt.show()


