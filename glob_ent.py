'''
This module studies radiation at a global scale.
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import calendar as cal
from ent_datasets import months_in_year

# Load datasets
grid = np.load('datasets/earth_grid.npy')
lat, lon = grid[:,0,0], grid[0,:,1]

wvlen = np.load('datasets/wvlen.npy')           # nm
wvlen_lres = np.load('datasets/wvlen_lres.npy') # nm
wvnum = np.load('datasets/wvnum.npy')           # cm^-1
wvlen_num = 1.0e7/wvnum                         # nm

def loadflux(month, re='r'):
    '''
    Loads SW and LW radiation flux datasets which correspond to given month.
    
    SW radiation is corrected by a factor of 1000 if not lres.
    SW radiation units: uW cm^-2 sr^-1 nm^-1
    LW radiation units: W cm^-2 sr^-1 cm
    
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    Parameter re takes the values 'r' and 'e' for radiation and entropy 
    fluxes respectively.
    '''
    swl = np.load('datasets/flux/lres_sw{0}{1}.npy'.format(re,month))
    lw = np.load('datasets/flux/lw{0}{1}.npy'.format(re,month))
    lwc = np.load('datasets/flux/clr_lw{0}{1}.npy'.format(re,month))
    return np.array([swl, lw, lwc])

titles = ['SW (low res)', 'LW', 'Clear sky LW']
rad = {'swl':0, 'lw':1, 'lwc':2}

# Plot functions
def plot_map(xx, yy, data, u):
    '''
    Plots a global map of the data provided.
    
    Parameters correspond to londitude and latitude, data to be plotted and 
    units.
    '''
    m = Basemap()
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,data)
    m.colorbar(location='bottom', label='Flux '+u, pad="10%")
    parallels = np.arange(-90, 90, 22.5)
    meridians = np.arange(-180,180, 30)
    m.drawparallels(parallels, labels=[1,0,0,0])
    m.drawmeridians(meridians, labels=[0,0,0,1])

def plot_flux(month, flux, re='r', info=False):
    '''
    Plots global maps of given SW and LW fluxes.
    
    Parameter flux is an array of fluxes.
    '''
    xx, yy, shflux = shift_grid(flux)
    if info: return xx, yy, shflux
    
    s, u = is_re(re)
    m = calendar(month)
    fig = plt.figure()
    fig.suptitle('Global {0} flux for {1} {2}'
                 .format(s, m[1], m[0]), fontsize=16)
    for i in range(len(shflux)):
        ax = fig.add_subplot(2,2,i+2)
        ax.set_title('{0} flux'.format(titles[i]), fontsize=14)
        plot_map(xx, yy, shflux[i], u)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.925)

def plot_diff(month1, month2, flux1, re='r', info=False):
    '''
    Plots the difference of fluxes between two months.
    If info=True, the difference is returned and not plotted.
    
    Parameter month2 corresponds to the month, month1 (with flux1) is compared 
    with.
    '''
    flux2 = loadflux(month2, re)
    fdiff = flux1-flux2
    xx, yy, shfdiff = shift_grid(fdiff)
    if info: return xx, yy, shfdiff
    
    s, u = is_re(re)
    m1, m2 = calendar(month1), calendar(month2)
    fig = plt.figure()
    fig.suptitle('Difference of {0} flux between {1} {2} and {3} {4}'
                 .format(s, m1[1], m1[0], m2[1], m2[0]), fontsize=16)
    
    for i in range(len(fdiff)):
        ax = fig.add_subplot(2,2,i+2)
        ax.set_title('{0} flux difference'.format(titles[i]), fontsize=14)
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
    Returns useful information about array arr (eg. extreme values).
    '''
    return {'shape':arr.shape, 'max':np.nanmax(arr), 'min':np.nanmin(arr),
            'nans':np.where(np.isnan(arr))[0].size}

def calendar(month):
    '''
    Converts a string of the form yymm (y for year, m for month) to yyyy and 
    the name of the corresponding month m (eg. input '0001' returns 
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
    Returns 'radiation' and 'entropy' to be used in titles and plot labels for 
    re = 'r' and 'e' respectively, with the corresponding flux units.
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
    Plots SW & LW radiation and entropy fluxes across the globe. It can also 
    call function sumup for the flux data of the month.
    
    Parameter diff is a list of months to be compared with month, gmap and 
    info plot the relevant maps and provide information on the data 
    respectively, if True.
    
    For all types of radiation, info calls sumup(), then returns the flux data 
    for given month and re, then returns  the difference in flux data from 
    months in the diff list.
    '''
    flux = loadflux(month, re)
    if gmap:
        plot_flux(month, flux, re, info=False)
        if diff is not None:
            for m in diff:
                plot_diff(month, m, flux, re, info=False)
    
    if info:
        suma, dflux, ddiff = {}, [], []
        for f in range(len(flux)):
            suma[rad.keys()[rad.values().index(f)]] = sumup(flux[f])
        dflux.append(plot_flux(month, flux, re, info=True)[2])
        for m in diff:
            ddiff.append(plot_diff(month, m, flux, re, info=True)[2])
        return suma, dflux, ddiff

def spec_analysis(month, re='r', radtype='swl', gmap=True, info=False):
    '''
    Plots monthly comparisons of specific type of flux across the globe 
    (SW, LW, Clear sky LW).
    
    Parameter radtype ('swl', 'lw', 'lwc') determines the type of flux. 
    Parameter info returns the latitude and longitude grid and the differences 
    between flux data for neighbouring months.
    '''
    flux = np.array([loadflux(month, re)[rad[radtype]]])
    s, u = is_re(re)
    neighbours = neighbour_months(month)
    m = calendar(month)
    
    fig = plt.figure()
    fig.suptitle('{0} {1} flux comparison for {2} {3}'
                 .format(titles[rad[radtype]], s, m[1], m[0]), fontsize=16)
    
    ax = fig.add_subplot(2,2,1)
    ax.set_title('Flux in {0} {1}'.format(m[1], m[0]), fontsize=14)
    xx, yy, shflux = shift_grid(flux)
    if gmap: plot_map(xx, yy, shflux[0], u)
    
    fdiffs = [(xx, yy)]
    for n in range(len(neighbours)):
        ax = fig.add_subplot(2,2,n+2)
        m2 = calendar(neighbours[n])
        ax.set_title('Difference with {0} {1}'.format(m2[1], m2[0]), fontsize=14)
        flux2 = np.array([loadflux(neighbours[n], re)[rad[radtype]]])
        fdiff = flux - flux2
        xx, yy, fdiff = shift_grid(fdiff)
        fdiffs.append(fdiff[0])
        if gmap: plot_map(xx, yy, fdiff[0], u)
    plt.tight_layout()
    #plt.subplots_adjust(top=0.925)
    if info: return fdiffs

def analyse_month_mini(month, re='r', radtype='swl', export=False):
    '''
    Plots maps of given month for both 2000 and 2099 as well as a map with 
    their difference.
    '''
    if month[:2]=='00':
        month2, yr = '99'+month[2:], ['2000', '2099']
    elif month[:2]=='99': month2, yr = '00'+month[2:], ['2099', '2000']
    flux = np.array([loadflux(month, re)[rad[radtype]]])
    flux2 = np.array([loadflux(month2, re)[rad[radtype]]])
    s, u = is_re(re)
    m = calendar(month)
    
    fig = plt.figure()
    fig.suptitle('{0} {1} flux comparison for {2}'
                 .format(titles[rad[radtype]], s, m[1]), fontsize=16)
    
    ax = fig.add_subplot(2,2,1)
    ax.set_title('Flux in {0} {1}'.format(m[1], yr[0]), fontsize=14)
    xx, yy, shflux = shift_grid(flux)
    plot_map(xx, yy, shflux[0], u)
    
    ax2 = fig.add_subplot(2,2,2)
    ax2.set_title('Flux in {0} {1}'.format(m[1], yr[1]), fontsize=14)
    xx, yy, shflux2 = shift_grid(flux2)
    plot_map(xx, yy, shflux2[0], u)
    
    axd = fig.add_subplot(2,2,3)
    axd.set_title('Difference in flux', fontsize=14)
    xx, yy, shfluxd = shift_grid(flux2-flux)
    plot_map(xx, yy, shfluxd[0], u)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.925)
    
    if export:
        fig.savefig('datasets/maps/{0}_{1}_{2}.png'
                    .format(re,radtype,month[2:]))
        plt.close(fig)

def export(year):
    '''
    Exports all maps from analyse_month_mini with given year as a parameter.
    '''
    for m in months_in_year(2000):
        analyse_month_mini(m, re='r', radtype='swl', export=True)
        analyse_month_mini(m, re='r', radtype='lwc', export=True)
        analyse_month_mini(m, re='r', radtype='lw', export=True)
        analyse_month_mini(m, re='e', radtype='swl', export=True)
        analyse_month_mini(m, re='e', radtype='lwc', export=True)
        analyse_month_mini(m, re='e', radtype='lw', export=True)
        print 'DONE', m

#analyse_month_mini('0010', re='e', radtype='lw', export=False)
#plot_diff('9907','9908', loadflux('9907', 'e'), re='e')
suma, dflux, ddiff = analyse_month('9902', re='r', diff=['9902','0004'], gmap=False, info=True)
print suma
print dflux
print ddiff
#fdiffs = spec_analysis('9902', 'r', 'lw',gmap=True,info=True)
#print fdiffs
#m=months_in_year(2099)
#spec_analysis(m[8], 'e', 'lwc')
plt.show()


