'''
This module reads libradtran output and calculates material entropy.
'''

import numpy as np
import collections
import matplotlib.pylab as plt
import scipy.interpolate as sci
from glob_ent import find_nans

# Physical constants
c=2.998e8       # m s^-1
kb=1.38e-23     # kg m^2 s^-2 K^-1
h=6.626e-34     # J s = kg m^2 s^-1
sb=5.670e-8     # W m^-2 K^-4 = kg s^-3 K^-4
cp=1005         # J kg^-1 K^-1 = m^2 s^-2 K^-1 (specific cp of air)
sid=60*60*24.0  # sec (seconds in a day)
pi = np.pi      # sr (solid angle assumed for isotropically diffused radiation)
omega0=6.87e-5  # sr (solid angle subtended by the Sun)

# Load data
def load_output(filename, radtype='sw', only_output=False):
    '''
    Reads libradtran output file. The physical quantities are:
    
    - wvl/wvn (wavelength, nm converts to wavenumber, m^-1)
    - z (altitude, km converts to m)
    - T (temperature, K)
    - rho_air (air density, kg m^-3)
    - edir (direct flux, W m^-2 nm^-1 converts to W m^-2 m) !!! kurudz has mW??
    - edn (upwards diffuse flux, same as edir)
    - eup (downwards diffuse flux, same as edir)
    - heat (heating rate, K day^-1 nm^-1 converts to K day^-1 m)
    '''
    fdir = "libradtran/{0}.out".format(filename)
    output = np.loadtxt(fdir, comments='#')
    znum = np.count_nonzero(output[:,0]==output[:,0][0])
    
    wvl = output[:,0][::znum]
    wvn = 1.0e9/wvl
    z = output[:,1][:znum]*1000
    T = output[:,2][:znum]
    rho_air = output[:,3][:znum]
    
    output = output[:,4:].reshape(len(wvl), znum, 4)
    units = 1.0e-9*wvl[:,None]*wvl[:,None]
    edir = output[:,:,0] * units
    edn = output[:,:,1] * units
    eup = output[:,:,2] * units
    heat = output[:,:,3] * units
    
    if only_output: return edir, edn, eup, heat 
    return wvn, z, T, rho_air, edir, edn, eup, heat

# Math and physics conversion functions
def ent_flux(wvn, rad, radtype='sw', angle='dir'):
    '''
    Converts intensity to entropy flux with units W m^-2 m K^-1.
    
    Taylor approximation is used as in ent_datasets.py.
    '''
    if angle=='dir': ang = omega0
    elif angle=='diff': ang = pi
    intens = rad/ang
    wvn = wvn[:,None]
    
    y = intens/(2*h*c*c*wvn**3)
    np.seterr(all='ignore')
    ent1 = np.where(y>=0.01, (1+y)*np.log(1+y), (1+y)*(y - y*y/2 + y**3/3))
    # Since y > -1e10, all y < 0 are converted to 0 in ent2
    ent2 = np.where(y>0.0, -y*np.log(y), 0.0) 
    ent = np.where(y!=0.0, 2*kb*c*wvn*wvn*(ent1+ent2)*ang, 0.0)
    np.seterr(all='warn')
    return ent

def deriv(z, rad):
    '''
    Evaluates derivative of radiation flux or entropy with respect to altitude 
    using B-splines (1D z array and 2D rad array).
    '''
    derv=np.zeros(rad.shape)
    for i in range(rad.shape[0]):
        tck = sci.splrep(z, rad[i,:])
        derv[i,:] = sci.splev(z, tck, der=1)
    return derv

def spect_int(wvn, rad):
    '''
    Evaluates the spectral integral of given radiation flux or entropy over 
    given wavenumber range (1D wvn array and 2D rad array).
    '''
    intg = np.zeros(rad.shape[1])
    for i in range(rad.shape[1]):
        intg[i] = -np.trapz(rad[:,i], wvn)
    return intg

# Analysis
def heat_check(z, rho_air, edir, edn, eup):
    '''
    Estimates the libradtran heat input in units K day^-1 m for comparison.
    '''
    const = sid/cp/rho_air[None,:]
    check = (deriv(z, edir)+deriv(z, edn)+deriv(z, eup))*const
    return check

def sdot(quants, radtype='sw', debug=False):
    '''
    Returns material, radiation and total entropy rates and a heat check.
    Entropy rates are in units W m^-3 K^-1 m.
    
    Parameter quants is the output of the load_output function when 
    only_output=False.
    '''
    wvn, z, T, rho_air, edir, edn, eup, heat = quants
    ent_dir = ent_flux(wvn, edir, radtype=radtype, angle='dir')
    ent_dn = ent_flux(wvn, edn, radtype=radtype, angle='diff')
    ent_up = ent_flux(wvn, eup, radtype=radtype, angle='diff')
    
    sdotmat = cp/sid * rho_air[None,:]/T[None,:] * heat
    sdotrad = deriv(z, ent_up) - deriv(z, ent_dn) - deriv(z, ent_dir)
    sdot = sdotmat + sdotrad
    check = heat_check(z, rho_air, edir, edn, eup)
    
    if debug:
        print ent_dn.size,'\n'
        print 'ent_dir', len(ent_dir[ent_dir==0])#, ' / ', len(edir[edir<=0])
        print '--------------------'
        print 'ent_dn', len(ent_dn[ent_dn==0])
        print '--------------------'
        print 'ent_up', len(ent_up[ent_up==0])
        print '--------------------'
        print 'sdotrad', len(sdotrad[sdotrad<0])
        print sdotrad
    
    return sdotmat, sdotrad, sdot, check

# Export output

def flux_output(filename, radtype='sw'):
    '''
    Exports flux output into a .txt file.

    F - radiation flux integrated over wavelength
    J - entropy flux integrated over wavelength
    Q - heat
    sdot - entropy rate
    '''
    quants = load_output(filename, radtype=radtype, only_output=False)
    wvn, z, T, rho_air, edir, edn, eup, heat = quants
    sdots = sdot(quants, radtype=radtype, debug=False)
    F_dir = spect_int(wvn, edir)
    F_up = spect_int(wvn, eup)
    F_dn = spect_int(wvn, edn)
    J_dir = spect_int(wvn, ent_flux(wvn, edir, radtype=radtype, angle='dir'))
    J_up = spect_int(wvn, ent_flux(wvn, eup, radtype=radtype, angle='diff'))
    J_dn = spect_int(wvn, ent_flux(wvn, edn, radtype=radtype, angle='diff'))
    s = spect_int(wvn, sdots[2])
    s_r = spect_int(wvn, sdots[1])
    s_m = spect_int(wvn, sdots[0])
    Q = spect_int(wvn, heat)
    Q_check = spect_int(wvn, sdots[3])
    
    out = open("libradtran/entropy_budget_{0}.txt".format(radtype), 'w')
    out.write('{0:>10} {1:>11} {2:>11} {3:>11} {4:>11} {5:>11} {6:>11} {7:>11} {8:>11} {9:>11} {10:>11} {11:>11}\n'.format('z ', 'F_dir ', 'F_up ', 'F_dn ', 'J_dir ', 'J_up ', 'J_dn ', 'sdot ', 'sdot_r ', 'sdot_m ', 'Q ', 'Q_check '))
    out.write('{0:>10} {1:>11} {2:>11} {3:>11} {4:>11} {5:>11} {6:>11} {7:>11} {8:>11} {9:>11} {10:>11} {11:>11}\n'.format('(km)', '(W m-2)', '(W m-2)', '(W m-2)', '(W m-2 K-1)', '(W m-2 K-1)', '(W m-2 K-1)', '(W m-3 K-1)', '(W m-3 K-1)', '(W m-3 K-1)', '(K day-1)', '(K day-1)'))
    for i in range(len(z)):
        out.write('{0:>10.3f} {1:>11.2E} {2:>11.2E} {3:>11.2E} {4:>11.2E} {5:>11.2E} {6:>11.2E} {7:>11.2E} {8:>11.2E} {9:>11.2E} {10:>11.2E} {11:>11.2E}\n'.format(z[i], F_dir[i], F_up[i], F_dn[i], J_dir[i], J_up[i], J_dn[i], s[i], s_r[i], s_m[i], Q[i], Q_check[i]))
    
    rad = J_up[0]-J_dn[0]-J_dir[0]
    mat = -(F_up[0]-F_dn[0]-F_dir[0])/T[0]
    net = rad + mat
    out.write('\nSurface entropy (W m-2 K-1): rad = {0}, mat = {1}, net = {2}'
              .format(rad, mat, net))
    
    out.close()

flux_output('0solar_rep', radtype='sw')







