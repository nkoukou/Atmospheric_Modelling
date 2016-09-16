'''
This module reads libradtran output and calculates material entropy.
'''

import numpy as np
import collections
import matplotlib.pylab as plt
import scipy.interpolate as sci

# Physical constants
c=2.998e8       # m s^-1
kb=1.38e-23     # kg m^2 s^-2 K^-1
h=6.626e-34     # J s = kg m^2 s^-1
sb=5.670e-8     # W m^-2 K^-4 = kg s^-3 K^-4
cp=1005         # J kg^-1 K^-1 = m^2 s^-2 K^-1 (specific cp of air)
sid=60*60*24.0  # sec (seconds in a day)
pi = np.pi      # sr (solid angle assumed for isotropically diffused radiation)
omega0=6.87e-5  # sr (solid angle subtended by the Sun)

def load_output(filename, radtype='sw', only_output=False):
    '''
    Reads libradtran output file. The physical quantities are:
    
    - wvl (wavelength, nm)
    - z (altitude, km converts to m)
    - T (temperature, K)
    - rho_air (air density, kg m^-3)
    - output (3D array containing the variable quantities listed below)
        - edir (direct intensity)
        - edn (upwards diffuse intensity)
        - eup (downwards diffuse intensity)
        - heat (heating rate, K day^-1 nm^-1)
        All intensities are read in units mW m^-2 nm^-1 which convert into 
        W m^-2 m. heat is transformed to K day^-1 m.
    '''
    fdir = "libradtran/{0}.out".format(filename)
    output = np.loadtxt(fdir, comments='#')
    znum = np.count_nonzero(output[:,0]==output[:,0][0])
    
    wvl = output[:,0][::znum]
    z = output[:,1][:znum]*1000
    T = output[:,2][:znum]
    rho_air = output[:,3][:znum]
    
    output = output[:,4:].reshape(len(wvl), znum, 4)
    units = 1.0e-12*wvl[:,None]*wvl[:,None]
    edir = output[:,:,0] * units
    edn = output[:,:,1] * units
    eup = output[:,:,2] * units
    heat = output[:,:,3] * units * 1000
    
    if only_output: return edir, edn, eup, heat 
    return wvl, z, T, rho_air, edir, edn, eup, heat

def ent_flux(wvl, rad, radtype='sw', angle='dir'):
    '''
    Transforms intensity to entropy flux with units W m^-2 m K^-1.
    '''
    if angle=='dir': ang = omega0
    elif angle=='diff': ang = pi
    intens = rad/ang
    wvn = 1.0e9/wvl[:,None]
    
    y = intens/(2*h*c*c*wvn**3)
    np.seterr(all='ignore')
    ent1 = np.where(y>=0.01, (1+y)*np.log(1+y), (1+y)*(y - y*y/2 + y**3/3))
    ent2 = np.where(y!=0.0, -y*np.log(y), 0.0)
    ent = np.where(y!=0.0, 2*kb*c*wvn*wvn*(ent1+ent2)*ang, 0.0)
    np.seterr(all='warn')
    return ent

def interpol(z, rad):
    '''
    Performs Lagrange interpolation on 1D z array and 2D radiation array
    '''
    inter=np.zeros(rad.shape)
    for i in range(len(rad)):
        tck = sci.splrep(z, rad[i,:])
        inter[i,:] = sci.splev(z, tck, der=1)
    return inter

def heat_check(z, edir, edn, eup):
    '''
    Estimates the libradtran heat input in units K day^-1 m for comparison.
    '''
    const = sid/cp/rho_air[None,:]
    check = (interpol(z, edir)+interpol(z, edn)+interpol(z, eup))*const
    return check

#wvl, z, T, rho_air, edir, edn, eup, heat = load_output('testSO')
#inter = interpol(z, edir)
#plt.plot(z, edir[0, :], z, inter[0, :])
#plt.show()







