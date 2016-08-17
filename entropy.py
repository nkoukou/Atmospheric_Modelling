'''
This module transfers 'entropy.pro' into python code.
'''

import pickle as pk
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm

# Load datasets
dataset = open("datasets/gm_dict.pkl",'rb')
gm_dict = pk.load(dataset)
dataset.close()

wvlen = gm_dict['wvl']  # nm
wvnum = gm_dict['wvn']  # cm^-1
wvlen_num = 1.0e7/wvnum # nm

sw_gm00_rad = gm_dict['sw00'] # !!!
sw_gm99_rad = gm_dict['sw99'] #
lw_gm00_rad = gm_dict['lw00'] #
lw_gm99_rad = gm_dict['lw99'] #


# Physical constants
c=2.998e8   # m s^-1
kb=1.38e-23 # kg m^2 s^-2 K^-1
h=6.626e-34 # J s = kg m^2 s^-1

a = [[1,2,3],[2,3,4]]



# Choose year 2000 or 2099
def choose_year(year):
    '''
    Determines all variables depending on the year of choice.
    Involves correction for sw rad by a factor of 1000.
    '''
    if year==2000:
        sw_gm00_tavg, lw_gm00_tavg = 0, 0
        for month in range(12):
            sw_gm00_tavg += sw_gm00_rad[:,month]
            lw_gm00_tavg += lw_gm00_rad[:,month]
        sw_gm00_tavg, lw_gm00_tavg = sw_gm00_tavg/12.0, lw_gm00_tavg/12.0
        return sw_gm00_rad, lw_gm00_rad, sw_gm00_tavg, lw_gm00_tavg
    elif year==2099:
        sw_gm99_tavg, lw_gm99_tavg = 0, 0
        for month in range(12):
            sw_gm99_tavg += sw_gm99_rad[:,month]
            lw_gm99_tavg += lw_gm99_rad[:,month]
        sw_gm99_tavg, lw_gm99_tavg = sw_gm99_tavg/12.0, lw_gm99_tavg/12.0
        return sw_gm99_rad, lw_gm99_rad, sw_gm99_tavg, lw_gm99_tavg

# Conversion functions
# Fix units in docstring !!!
def radtoent(rad,radtype='sw'):
    '''
    Converts array of radiation intensity to array of entropy
    '''
    if radtype=='sw':
        wvl=wvlen
        iconst=1.0e-11*wvl*wvl
    elif radtype=='lw':
        wvl=wvlen_num
        iconst=1.0e2
    wvn = 1.0e9/wvl
    econst=1.0e-6*wvn*wvn
    
    intens = iconst*rad
    y = intens/(2*h*c*c*wvn**3)
    ent = econst*2*kb*c*wvn*wvn*((1+y)*np.log(1+y)-y*np.log(y))
    return ent

def radtorad(rad,radtype='sw'):
    '''
    Converts radiation to units W m-2 Sr-1 nm-1
    '''
    if radtype=='sw':
        wvln=wvlen
        rconst= 1.0e-2
    elif radtype=='lw':
        wvln=wvnum
        rconst= 1.0e-3*wvnum*wvnum
    
    return rconst*rad

def rad_flux(month,year=2000,radtype='sw'):
    '''
    Integrates radiation over wavelength to get radiative flux
    '''
    if radtype=='sw':
        wvl = wvlen
        rad = choose_year(year)[0][:,month]
    elif radtype=='lw':
        wvl = wvlen_num
        rad = choose_year(year)[1][:,month]
    n = len(wvl)-1
    rad = radtorad(rad,radtype)
    
    flux = rad[0]*(wvl[0]-wvl[1])+rad[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*rad[i]*(wvl[i-1]-wvl[i+1])
    return flux

def ent_flux(month,year=2000,radtype='sw'):
    '''
    Integrates entropy over wavelength to get radiative flux
    '''
    if radtype=='sw':
        wvl = wvlen
        rad = choose_year(year)[0][:,month]
    elif radtype=='lw':
        wvl = wvlen_num
        rad = choose_year(year)[1][:,month]
    n = len(wvl)-1
    ent = radtoent(rad,radtype)
    
    flux = ent[0]*(wvl[0]-wvl[1])+ent[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*ent[i]*(wvl[i-1]-wvl[i+1])
    return flux

# Plot functions
# Redesign this, tex font bigger, units for yaxis, docstrings,
# 1e-7 overlaps with title, are graphs the same as in IDL?
# plot LW on top of SW in ent_rad_diff !!!
def plot_year_rad(year):
    fig, ((ax_swrm,ax_swrd),(ax_lwrm,ax_lwrd))\
    = plt.subplots(2,2,figsize=(9,7))
    
    sw_gm_rad, lw_gm_rad, sw_gm_tavg, lw_gm_tavg = choose_year(year)
    
    colors = cm.rainbow(np.linspace(0, 1, 12))
    
    for month in range(12):
        ax_swrm.plot(wvlen,sw_gm_rad[:,month],
                    '-',color=colors[month],lw=0.2)
        ax_swrd.plot(wvlen,sw_gm_rad[:,month]-sw_gm_tavg,
                    '-',color=colors[month],lw=0.2)
        ax_lwrm.plot(wvnum, lw_gm_rad[:,month],
                    '-',color=colors[month],lw=0.2)
        ax_lwrd.plot(wvnum, lw_gm_rad[:,month]-lw_gm_tavg,
                    '-',color=colors[month],lw=0.2)
    
    fig.suptitle('Global mean radiation for all months in year '
                +str(year), size=16)
    ax_swrm.set_title('Mean SW rad',size=12)
    ax_swrm.set_xlabel('Wavelength ($nm$)')
    ax_swrd.set_title('Deviation of SW rad from annual mean',size=12)
    ax_swrd.set_xlabel('Wavelength ($nm$)')
    ax_lwrm.set_title('Mean LW rad',size=12)
    ax_lwrm.set_xlabel('Wavenumber ($cm^{-1}$)')
    ax_lwrd.set_title('Deviation of LW rad from annual mean',size=12)
    ax_lwrd.set_xlabel('Wavenumber ($cm^{-1}$)')
    plt.tight_layout()

def plot_month_rad(month):
    fig, ((ax_swr_ann,ax_swrdiff),(ax_lwr_ann,ax_lwrdiff))\
    = plt.subplots(2,2,figsize=(9,7))
    ax_swr_ann.plot(wvlen,sw_gm00_rad[:,month],
                '-',lw=0.2,color='b', label='2000')
    ax_swr_ann.plot(wvlen,sw_gm99_rad[:,month],
                '-',lw=0.2,color='g', label='2099')
    ax_swrdiff.plot(wvlen,
      sw_gm99_rad[:,month]-sw_gm00_rad[:,month],
      '-',lw=0.2,color='r')
    
    ax_lwr_ann.plot(wvnum, lw_gm00_rad[:,month],
                '-',lw=0.2,color='b', label='2000')
    ax_lwr_ann.plot(wvnum, lw_gm99_rad[:,month],
                '-',lw=0.2,color='g', label='2099')
    ax_lwrdiff.plot(wvnum,
      lw_gm99_rad[:,month]-lw_gm00_rad[:,month],
      '-',lw=0.2,color='r')
    
    fig.suptitle('Global mean radiation for month '+str(month+1), size=16)
    ax_swr_ann.set_title('SW rad for years 2000 and 2099',size=12)
    ax_swr_ann.set_xlabel('Wavelength ($nm$)')
    sw_leg = ax_swr_ann.legend(loc='upper right')
    for legobj in sw_leg.legendHandles:
        legobj.set_linewidth(2.0)
    ax_swrdiff.set_title('Difference of SW rad between 2000 and\
                          2099',size=12)
    ax_swrdiff.set_xlabel('Wavelength ($nm$)')
    ax_lwr_ann.set_title('LW rad for years 2000 and 2099',size=12)
    ax_lwr_ann.set_xlabel('Wavenumber ($cm^{-1}$)')
    lw_leg = ax_lwr_ann.legend(loc='upper right')
    for legobj in lw_leg.legendHandles:
        legobj.set_linewidth(2.0)
    ax_lwrdiff.set_title('Difference of LW rad between 2000 and\
                          2099',size=12)
    ax_lwrdiff.set_xlabel('Wavenumber ($cm^{-1}$)')
    plt.tight_layout()

def plot_ent_month(month):
    fig, ((ax_sw00,ax_lw00),(ax_sw99,ax_lw99))\
    = plt.subplots(2,2,figsize=(9,7))
    ent_sw00 = radtoent(sw_gm00_rad[:,month],'sw')
    ent_sw99 = radtoent(sw_gm99_rad[:,month],'sw')
    ent_lw00 = radtoent(lw_gm00_rad[:,month],'lw')
    ent_lw99 = radtoent(lw_gm99_rad[:,month],'lw')
    
    ax_sw00.plot(wvlen,ent_sw00)
    ax_sw99.plot(wvlen,ent_sw99)
    ax_lw00.plot(wvlen_num,ent_lw00)
    ax_lw99.plot(wvlen_num,ent_lw99)

    fig.suptitle('Global mean entropy for month '+str(month+1), size=16)
    ax_sw00.set_title('SW entropy for year 2000')
    ax_sw00.set_xlabel('Wavelength ($nm$)')
    ax_sw99.set_title('SW entropy for year 2099')
    ax_sw99.set_xlabel('Wavelength ($nm$)')
    ax_lw00.set_title('LW entropy for year 2000')
    ax_lw00.set_xlabel('Wavelength ($nm$)')
    ax_lw99.set_title('LW entropy for year 2099')
    ax_lw99.set_xlabel('Wavelength ($nm$)')
    plt.tight_layout()

def plot_ent_rad_diff(month,radtype='sw'):
    if radtype=='sw':
        wvl = wvlen
        rad00 = choose_year(2000)[0][:,month]
        rad99 = choose_year(2099)[0][:,month]
        s='SW'
    elif radtype=='lw':
        wvl = wvlen_num
        rad00 = choose_year(2000)[1][:,month]
        rad99 = choose_year(2099)[1][:,month]
        s='LW'
    wvlog = np.log(wvl)
    rad00, rad99 = radtorad(rad00,radtype), radtorad(rad99,radtype)
    ent00, ent99 = radtoent(rad00,radtype), radtoent(rad99,radtype)
    
    fig, ((ax_rann,ax_eann),(ax_rdiff,ax_ediff))\
    = plt.subplots(2,2,figsize=(9,7))
    ax_rann.plot(wvlog,wvl*rad00,'-',lw=0.2,color='b', label='2000')
    ax_rann.plot(wvlog,wvl*rad99,'-',lw=0.2,color='g', label='2099')
    ax_eann.plot(wvlog,wvl*ent00,'-',lw=0.2,color='b', label='2000')
    ax_eann.plot(wvlog,wvl*ent99,'-',lw=0.2,color='g', label='2099')
    ax_rdiff.plot(wvlog,wvl*(rad99-rad00))
    ax_ediff.plot(wvlog,wvl*(ent99-ent00))
    
    fig.suptitle('Entropy and radiation difference between 2000 and 2099',
                size=16)
    ax_rann.set_title(s+' radiation',size=12)
    rleg = ax_rann.legend(loc='upper right')
    for legobj in rleg.legendHandles:
        legobj.set_linewidth(2.0)
    ax_eann.set_title('Entropy of '+s+' radiation',size=12)
    eleg = ax_eann.legend(loc='upper right')
    for legobj in eleg.legendHandles:
        legobj.set_linewidth(2.0)
    ax_rdiff.set_title('Difference in radiation between 2000 and 2099',
                      size=12)
    ax_ediff.set_title('Difference in entropy between 2000 and 2099',size=12)
    ax_rann.set_xlabel('$log\ \lambda$')
    ax_eann.set_xlabel('$log\ \lambda$')
    ax_rdiff.set_xlabel('$log\ \lambda$')
    ax_ediff.set_xlabel('$log\ \lambda$')
    plt.tight_layout()

def plot_swvslw(month,year):
    srad, lrad = choose_year(year)[0][:,month], choose_year(year)[1][:,month]
    
    fig, ax = plt.subplots(1,1,figsize=(9,7))
    ax.plot(np.log(wvlen),wvlen*srad,'-',lw=0.5,color='b', label='SW')
    ax.plot(np.log(wvlen_num),wvlen_num*lrad,'-',lw=0.5,color='r',\
           label='LW')
    ax.set_title('SW and LW radiation for year '+str(year))
    ax.set_xlabel('$log\ \lambda$')
    ax.legend(loc='upper right')

plot_year_rad(2000)
plot_ent_rad_diff(9,radtype='sw')
plt.show()



