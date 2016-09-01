'''
This module transfers 'entropy.pro' into python code. It analyses global mean 
data provided by OSSE.

I - radiation intensity; F - I integrated over wavelength
L - entropy intensity; J - L integrated over wavelength
'''

import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker as mtk
import matplotlib.cm as cm
import calendar as cal

# Load datasets
# SW rad is corrected by a factor of 1000
wvlen = np.load('datasets/wvlen.npy') # nm
wvnum = np.load('datasets/wvnum.npy') # cm^-1
wvlen_num = 1.0e7/wvnum               # nm

def loadgmrad(year):
    '''
    Loads global mean SW and LW radiation intensity datasets which correspond 
    to given year.
    
    SW radiation is corrected by a factor of 1000 if not lres.
    SW radiation units: uW cm^-2 sr^-1 nm^-1
    LW and clear sky (clr) LW radiation units: W cm^-2 sr^-1 cm
    '''
    if year==2000: yr='00'
    elif year==2099: yr='99'
    gm_sw_rad = 1000*np.load('datasets/gm_sw'+yr+'.npy')
    gm_lw_rad = np.load('datasets/gm_lw'+yr+'.npy')
    gm_clr_lw_rad = np.load('datasets/gm_clr_lw'+yr+'.npy')
    return gm_sw_rad, gm_lw_rad, gm_clr_lw_rad

titles = ['SW', 'LW', 'clear sky LW']

# Physical constants (accuracy to 4sf)
c=2.998e8    # m s^-1
kb=1.381e-23 # kg m^2 s^-2 K^-1
h=6.626e-34  # J s = kg m^2 s^-1

# Utility functions
def tavg(rad_ent):
    '''
    Calculates annual average of radiation or entropy intensity.
    '''
    return rad_ent.mean(axis=1)

def choose_band(title, re='r'):
    '''
    Chooses wvlen/wvnum and units depending on the type of radiation intensity.
    '''
    if title=='SW':
        wvln, xu = wvlen, '$Wavelength\ (nm)$'
        if re=='r': yu = '$I\ (\mu W\ cm^{-2}\ sr^{-1}\ nm^{-1})$'
        elif re=='e': yu = '$L\ (\mu W\ cm^{-2}\ sr^{-1}\ nm^{-1}\ K^{-1})$'
    elif title=='LW' or title=='clear sky LW':
        wvln, xu = wvnum, '$Wavenumber\ (cm^{-1})$'
        if re=='r': yu = '$I\ (\mu W\ cm^{-2}\ sr^{-1}\ nm^{-1})$'
        elif re=='e': yu = '$L\ (\mu W\ cm^{-2}\ sr^{-1}\ nm^{-1}\ K^{-1})$'
    return wvln, xu, yu

def is_re(re):
    '''
    Returns 'radiation' and 'entropy' to be used in titles and plot labels for 
    re = 'r' and 'e' respectively.
    '''
    if re=='r': s = 'radiation'
    elif re=='e': s = 'entropy'
    return s

# Conversion functions
def radtorad(rad,radtype='SW'):
    '''
    Converts radiation to units W m^-2 sr^-1 nm^-1
    '''
    if radtype=='SW':
        rconst = 1.0e-2
    elif radtype=='LW' or radtype=='clear sky LW':
        rconst = 1.0e-3*wvnum*wvnum
        rconst = rconst[:,None]
    return rconst*rad

def radtoent(rad,radtype='SW'):
    '''
    Converts array of radiation intensity to array of entropy according to 
    formula 8 in 'Wu and Liu: Radiation Entropy Flux, published 14/05/2010'. 
    Approximates ln(1+y) ~ y - y*y/2 + y**3/3 (Maclaurin) for y < 1.0e-2 to 
    account for miscalculation of the np.log function (1.0e-2 is arbitrary).
    
    Needs radiation in W m^-2 sr^-1 m.
    Returns entropy in mW m^-2 sr^-1 K^-1 nm^-1.
    '''
    if radtype=='SW':
        wvl = wvlen
        iconst = 1.0e-11*wvl*wvl
        iconst = iconst[:,None]
    elif radtype=='LW' or radtype=='clear sky LW':
        wvl = wvlen_num
        iconst = 1.0e2
    wvn = 1.0e9/wvl
    econst = 1.0e-6*wvn*wvn
    econst, wvn = econst[:,None], wvn[:,None]
    
    intens = iconst*rad
    y = intens/(2*h*c*c*wvn**3)
    
    ent1 = np.where(y>=0.01, (1+y)*np.log(1+y), (1+y)*(y - y*y/2 + y**3/3))
    ent2 = -y*np.log(y)
    ent = econst*2*kb*c*wvn*wvn*(ent1+ent2)
    return ent

def rad_flux(year=2000,radtype='SW'):
    '''
    Integrates radiation intensity over wavelength to get radiative flux.
    The units are W m^-2 sr^-1.
    '''
    if radtype=='SW':
        wvl = wvlen
        rad = loadgmrad(year)[0]
    elif radtype=='LW':
        wvl = wvlen_num
        rad = loadgmrad(year)[1]
    elif radtype=='clear sky LW':
        wvl = wvlen_num
        rad = loadgmrad(year)[2]
    n = len(wvl)-1
    rad = radtorad(rad,radtype)
    
    flux = rad[0]*(wvl[0]-wvl[1])+rad[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*rad[i]*(wvl[i-1]-wvl[i+1])
    return flux

def ent_flux(year=2000,radtype='SW'):
    '''
    Integrates entropy intensity over wavelength to get radiative flux.
    The units are mW m^-2 sr^-1 K^-1.
    '''
    if radtype=='SW':
        wvl = wvlen
        rad = loadgmrad(year)[0]
    elif radtype=='LW':
        wvl = wvlen_num
        rad = loadgmrad(year)[1]
    elif radtype=='clear sky LW':
        wvl = wvlen_num
        rad = loadgmrad(year)[2]
    n = len(wvl)-1
    ent = radtoent(rad,radtype)
    
    flux = ent[0]*(wvl[0]-wvl[1])+ent[n]*(wvl[n-1]-wvl[n])
    for i in range(1,n-1):
        flux += 0.5*ent[i]*(wvl[i-1]-wvl[i+1])
    return flux

# Plot functions
def plot_year_rad(year):
    '''
    Plots radiation intensity against wavelength/wavenumber, and the radiation 
    deviation from the mean for the given year.
    '''
    gm_rad = loadgmrad(year)
    gm_tavg = []
    for i in range(len(gm_rad)):
        gm_tavg.append(tavg(gm_rad[i]))
        
    fig = plt.figure()
    fig.suptitle('Global mean radiation for all months in year '+str(year),
                 size=16)
    colors = cm.rainbow(np.linspace(0, 1, 12))
    for i in range(len(gm_rad)):
        axm = fig.add_subplot(3,2,2*i+1)
        axd = fig.add_subplot(3,2,2*i+2)
        title = titles[i]
        wvln, xu, yu = choose_band(title, re='r')
        for month in range(12):
            axm.plot(wvln, gm_rad[i][:,month], color=colors[month], lw=0.2)
            axd.plot(wvln, gm_rad[i][:,month]-gm_tavg[i],
                     color=colors[month], lw=0.2)
        axm.set_title('Mean '+title+' radiation', size=14)
        axd.set_title('Deviation of '+title+' radiation from annual mean',
                      size=14)
        axm.set_xlabel(xu)
        axd.set_xlabel(xu)
        axm.yaxis.set_major_formatter(mtk.FormatStrFormatter('%.2e'))
        axd.yaxis.set_major_formatter(mtk.FormatStrFormatter('%.2e'))
        axm.set_ylabel(yu)
        axd.set_ylabel(yu)
    plt.tight_layout(pad=-1.3, w_pad=-4.0, h_pad=-1.0)
    plt.subplots_adjust(top=0.925)

def plot_month(month, re='r'):
    '''
    Plots radiation or entropy intensities (re= 'r' or 'e') in 2000 and 2099 
    as well as the difference between the two years for given month.
    '''
    c = is_re(re)
    fig = plt.figure()
    fig.suptitle('Global mean '+c+' for '+cal.month_name[month],
                 size=16)
    for i in range(len(titles)):
        axm = fig.add_subplot(3,2,2*i+1)
        axd = fig.add_subplot(3,2,2*i+2)
        title = titles[i]
        if re=='r':
            gm00 = loadgmrad(2000)[i][:,month-1]
            gm99 = loadgmrad(2099)[i][:,month-1]
        elif re=='e':
            gm00 = radtoent(loadgmrad(2000)[i], title)[:,month-1]
            gm99 = radtoent(loadgmrad(2099)[i], title)[:,month-1]
        wvln, xu, yu = choose_band(title, re=re)
        axm.plot(wvln, gm00, color='b', lw=0.5, label='2000')
        axm.plot(wvln, gm99, color='g', lw=0.5, label='2099')
        axd.plot(wvln, gm99-gm00, color='r', lw=0.5)
        axm.set_title(title+' '+c+' for years 2000 and 2099', size=14)
        axd.set_title('Difference between 2000 and 2099', size=14)
        axm.set_xlabel(xu)
        axd.set_xlabel(xu)
        axm.yaxis.set_major_formatter(mtk.FormatStrFormatter('%.2e'))
        axd.yaxis.set_major_formatter(mtk.FormatStrFormatter('%.2e'))
        axm.set_ylabel(yu)
        axd.set_ylabel(yu)
        leg = axm.legend(loc='upper right')
        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)
    plt.tight_layout(pad=-1., w_pad=-6.2, h_pad=-1.0)
    plt.subplots_adjust(top=0.925)

def plot_flux(radtype='SW'):
    '''
    Plots radiation and entropy fluxes for 2000 and 2099 and compares the two 
    years, for given radiation type.
    '''
    ent00 = ent_flux(2000, radtype)
    ent99 = ent_flux(2099, radtype)
    rad00 = rad_flux(2000, radtype)
    rad99 = rad_flux(2099, radtype)
    
    fig, ((axr, axrd), (axe, axed)) = plt.subplots(2,2)
    fig.suptitle('Global mean radiation and entropy fluxes in 2000 and 2099',
                 size=16)
    months = range(len(rad00))
    
    axr.plot(months, rad00, color='b', label='2000')
    axr.plot(months, rad99, color='g', label='2099')
    axrd.plot(months, rad99-rad00, color='r')
    axe.plot(months, ent00, color='b', label='2000')
    axe.plot(months, ent99, color='g', label='2099')
    axed.plot(months, ent99-ent00, color='r')
    
    axr.set_title('Radiation flux', size=14)
    axrd.set_title('Radiation flux difference', size=14)
    axe.set_title('Entropy flux', size=14)
    axed.set_title('Entropy flux difference', size=14)
    axr.set_xlabel('Month (0 stands for January)')
    axrd.set_xlabel('Month (0 stands for January)')
    axe.set_xlabel('Month (0 stands for January)')
    axed.set_xlabel('Month (0 stands for January)')
    axr.set_ylabel('$F\ (W\ m^{-2}\ sr^{-1})$')
    axrd.set_ylabel('$\Delta F\ (W\ m^{-2}\ sr^{-1})$')
    axe.set_ylabel('$J\ (mW\ m^{-2}\ sr^{-1}\ K^{-1})$')
    axed.set_ylabel('$\Delta J\ (mW\ m^{-2}\ sr^{-1}\ K^{-1})$')
    leg = axr.legend(loc='upper right')
    leg = axe.legend(loc='upper right')
    plt.tight_layout()
    plt.subplots_adjust(top=0.925)

def plot_year_diff(month, radtype='SW'):
    '''
    Plots difference in wvl*rad or wvl*ent against log of wvl. The y axis is 
    not just radiation or entropy intensities so that the graph integrates 
    correctly over the logarithm of wavelength which corresponds to the x axis.
    '''
    
    rad00 = loadgmrad(2000)[titles.index(radtype)]
    rad99 = loadgmrad(2099)[titles.index(radtype)]
    ent00 = radtoent(rad00,radtype)[:,month-1]
    ent99 = radtoent(rad99,radtype)[:,month-1]
    rad00 = radtorad(rad00,radtype)[:,month-1]
    rad99 = radtorad(rad99,radtype)[:,month-1]
    
    if radtype=='SW':
        wvl = wvlen
    elif radtype=='LW' or radtype=='clear sky LW':
        wvl = wvlen_num
    wvlog = np.log(wvl)
    
    fig, ((axr,axrd),(axe,axed)) = plt.subplots(2,2)
    fig.suptitle(radtype+' entropy and radiation flux in '
                 +cal.month_name[month], size=16)
    axr.plot(wvlog,wvl*rad00,color='b', label='2000')
    axr.plot(wvlog,wvl*rad99,color='g', label='2099')
    axe.plot(wvlog,wvl*ent00,color='b', label='2000')
    axe.plot(wvlog,wvl*ent99,color='g', label='2099')
    axrd.plot(wvlog,wvl*(rad99-rad00),color='r')
    axed.plot(wvlog,wvl*(ent99-ent00),color='r')
    
    axr.set_title(radtype+' radiation',size=12)
    axe.set_title('Entropy of '+radtype+' radiation',size=12)
    axrd.set_title('Difference between 2000 and 2099', size=12)
    axed.set_title('Difference between 2000 and 2099',size=12)
    xu = '$log\ \lambda$'
    yur = '$\lambda I\ (W\ m^{-2}\ sr^{-1})$'
    yue = '$\lambda L\ (mW\ m^{-2}\ sr^{-1}\ K^{-1})$'
    axr.set_xlabel(xu); axr.set_ylabel(yur)
    axe.set_xlabel(xu); axe.set_ylabel(yue)
    axrd.set_xlabel(xu); axrd.set_ylabel(yur)
    axed.set_xlabel(xu); axed.set_ylabel(yue)
    axr.legend(loc='upper right')
    axe.legend(loc='upper right')
    plt.tight_layout()

def plot_type_diff(month, year):
    '''
    Plots wvl*rad for both SW and LW radiation against log of wvl for given 
    month and year in the same plot.
    '''
    srad = wvlen*radtorad(loadgmrad(year)[0],'SW')[:,month-1]
    lrad = wvlen_num*radtorad(loadgmrad(year)[1],'LW')[:,month-1]
    
    fig, ax = plt.subplots(1,1)
    ax.plot(np.log(wvlen),srad,color='b', label='SW')
    ax.plot(np.log(wvlen_num),lrad,color='g', label='LW')
    ax.set_title('SW and LW radiation flux for '+cal.month_name[month]+' '
                 +str(year))
    ax.set_xlabel('$log\ \lambda$')
    ax.set_ylabel('$F\ (W\ m^{-2}\ sr^{-1})$')
    ax.legend(loc='upper right')

#plot_swvslw(4,2099)
#plot_month(11,'e')
#plot_year_diff(12, 'LW')
#plot_type_diff(9,2099)
plt.show()



