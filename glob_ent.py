'''
This module studies radiation at a global scale.
'''

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import calendar as cal
import time

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
        rconst = 1.0e-2
    elif radtype=='lw':
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
    
    np.seterr(all='ignore')
    ent1 = np.where(y>=0.01, (1+y)*np.log(1+y), (1+y)*(y - y*y/2 + y**3/3))
    ent2 = -y*np.log(y)
    np.seterr(all='warn')
    ent = np.where(y!=0.0, econst*2*kb*c*wvn*wvn*(ent1+ent2), 0.0)
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
    s, u = is_re(re)
    
    xx, yy, sw_flux, lw_flux = shift_grid(sw_flux, lw_flux)
    
    fig = plt.figure(figsize=(12,10.5))
    fig.suptitle("Global "+s+" flux for "+month[1]+" "+month[0], fontsize=18)
    
    axsw = fig.add_subplot(211)
    axsw.set_title("SW flux", fontsize=15)
    m = Basemap()
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

def plot_diff(month1, month2, re='r', info=False):
    '''
    Parameter re takes the values 'r' and 'e' for radiation and entropy 
    fluxes respectively.
    
    This function plots the difference of SW and LW fluxes between two months.
    '''
    s, u = is_re(re)
    m1, m2 = calendar(month1), calendar(month2)
    
    fig = plt.figure(figsize=(12,10.5))
    fig.suptitle("Difference of "+s+" flux between "+m1[1]+" "
                 +m1[0]+" and "+m2[1]+" "+m2[0], fontsize=16)
    
    sw_flux1, lw_flux1 = flux_month(month1, re)
    sw_flux2, lw_flux2 = flux_month(month2, re)
    xx, yy, sw_flux1, lw_flux1 = shift_grid(sw_flux1, lw_flux1)
    sw_flux2, lw_flux2 = shift_grid(sw_flux2, lw_flux2)[2:]
    swf, lwf = sw_flux2-sw_flux1, lw_flux2-lw_flux1
    if info: return xx, yy, swf, lwf
    
    axsw = fig.add_subplot(211)
    axsw.set_title("SW flux difference", fontsize=14)
    m = Basemap()
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,swf)
    m.colorbar(location='bottom', label='Flux '+u)
    
    axlw = fig.add_subplot(212)
    axlw.set_title("LW flux difference", fontsize=14)
    m = Basemap()
    m.drawmapboundary()
    m.drawcoastlines()
    m.contourf(xx,yy,lwf)
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

def flux_month(month, re='r'):
    sw_rad, lw_rad = loadrad(month)
    if re=='r': s = rad_flux(sw_rad,'sw'); l = rad_flux(lw_rad,'lw')
    elif re=='e': s = ent_flux(sw_rad,'sw'); l = ent_flux(lw_rad,'lw')
    return s, l

def shift_grid(sw_flux, lw_flux):
    sw_flux, lons = shiftgrid(180.0, sw_flux, lon, start=False)
    lw_flux, lons = shiftgrid(180.0, lw_flux, lon, start=False)
    xx, yy = np.meshgrid(lons, lat)
    return xx, yy, sw_flux, lw_flux

# Main functions
def analyse_month(month, gmap=False, info=False):
    '''
    Parameter month must be of the form 'yymm' (e.g. '0001' corresponds
    to year 2000, month 01).
    
    This function plots SW & LW radiation and entropy fluxes across the 
    globe. It can also provide useful information about the datasets
    (e.g. extreme values)
    '''
    sw_rflux, lw_rflux = flux_month(month, re='r')
    sw_eflux, lw_eflux = flux_month(month, re='e')
    
    month = calendar(month)
    
    if gmap:
        plot_flux(month, sw_rflux, lw_rflux, 'r')
        plot_flux(month, sw_eflux, lw_eflux, 'e')
    
    if info:
        return sumup(sw_eflux), sumup(lw_eflux)

plot_diff('0009','0010')
#analyse_month('0004', True, True)
plt.show()


