'''
This module studies radiation at a global scale
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import calendar as cal

# Load datasets
grid = np.load('datasets/earth_grid.npy')
lat, lon = grid[:,0,0], grid[0,:,1]

wvlen = np.load('datasets/wvlen.npy') # nm
wvnum = np.load('datasets/wvnum.npy') # cm^-1
wvlen_num = 1.0e7/wvnum               # nm

def loadrad(month):
    '''
    Loads SW and LW radiation datasets which correspond to given month.
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    
    SW radiation is corrected by a factor of 1000.
    '''
    sw = 1000*np.load('datasets/sw%s.npy' % (month)) # uW cm^-2 sr^-1 nm^-1
    lw = np.load('datasets/lw%s.npy' % (month))      # W cm^-2 sr^-1 cm
    return sw, lw

# Physical constants
c=2.998e8   # m s^-1
kb=1.38e-23 # kg m^2 s^-2 K^-1
h=6.626e-34 # J s = kg m^2 s^-1

# Conversion functions
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
    ent = np.where(y!=0.0,
          econst*2*kb*c*wvn*wvn*((1+y)*np.log(1+y)-y*np.log(y)),0.0)
    return ent

def rad_flux(rad,radtype='sw'):
    '''
    Integrates radiation over wavelength to get radiative flux.
    The units are W m^-2 sr^-1.
    '''
    if radtype=='sw':
        wvl = wvlen
    elif radtype=='lw':
        wvl = wvlen_num
    n = len(wvl)-1
    wvl = wvl[:,None,None]
    rad = radtorad(rad,radtype)
    
    flux = rad[0]*(wvl[0]-wvl[1])+rad[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*rad[i]*(wvl[i-1]-wvl[i+1])
    return flux

def ent_flux(rad,radtype='sw'):
    '''
    Integrates entropy over wavelength to get radiative flux.
    The units are mW m^-2 sr^-1 K^-1.
    '''
    if radtype=='sw':
        wvl = wvlen
    elif radtype=='lw':
        wvl = wvlen_num
    n = len(wvl)-1
    wvl = wvl[:,None,None]
    ent = radtoent(rad,radtype)
    
    flux = ent[0]*(wvl[0]-wvl[1])+ent[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*ent[i]*(wvl[i-1]-wvl[i+1])
    return flux

# Plot functions
def plot_flux(month, sw_flux, lw_flux, re='r'):
    '''
    Parameter re takes the values 'r' and 'e' for radiation and entropy 
    fluxes respectively.
    
    This function plots global maps of given SW and LW fluxes.
    '''
    if re=='r': s, u = 'radiation', '($W\ m^{-2}\ sr^{-1}$)'
    elif re=='e': s, u = 'entropy', '($mW\ m^{-2}\ sr^{-1}\ K^{-1}$)'
    
    sw_flux, lons = shiftgrid(180.0, sw_flux, lon, start=False)
    lw_flux, lons = shiftgrid(180.0, lw_flux, lon, start=False)
    xx, yy = np.meshgrid(lons, lat)
    
    fig = plt.figure(figsize=(12,10.5))
    fig.suptitle("Global "+s+" flux for "+month[1]+" "+month[0], fontsize=18)
    
    axsw = fig.add_subplot(211)
    axsw.set_title("SW flux", fontsize=15)
    m = Basemap() #lat_0=0.0, lon_0=180.0
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,sw_flux)
    m.colorbar(location='bottom', label='Flux '+u)
    
    axlw = fig.add_subplot(212)
    axlw.set_title("LW flux", fontsize=15)
    m = Basemap()
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,lw_flux)
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

# Main functions
def analyse_month(month, info=False):
    '''
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    
    This function plots SW & LW radiation and entropy fluxes across the 
    globe. It can also provide useful information about the datasets
    (e.g. extreme values)
    '''
    sw_rad, lw_rad = loadrad(month)
    sw_rflux = rad_flux(sw_rad,'sw')
    sw_eflux = ent_flux(sw_rad,'sw')
    lw_rflux = rad_flux(lw_rad,'lw')
    lw_eflux = ent_flux(lw_rad,'lw')
    
    month = calendar(month)
    
    plot_flux(month, sw_rflux, lw_rflux, 'r')
    plot_flux(month, sw_eflux, lw_eflux, 'e')
    
    if info:
        return sumup(sw_eflux), sumup(lw_eflux)


analyse_month('0009')
plt.show()






