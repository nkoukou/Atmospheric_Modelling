'''
This module reads netcdf data files from the CCSM3 experiment and creates 
input for the libradtran code.
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pylab as plt

# Physical constants
kb=1.38e-23      # kg m^2 s^-2 K^-1
T0 = 273.15      # K
air0 = 1.2754e-3 # g cm^-3 (air density at T0 and P0=100kPa)

# Import netCDF files
m1 = "datasets/ccsm/control/b30.031.cam2.h0.0500-01.nc"
m2 = "datasets/ccsm/b30.031.cam2.h0.0359-01.nc"
m1 = Dataset(m1)
m2 = Dataset(m2)

# Testing data
def show_vars(filename):
    for var in filename.variables.keys():
        print var, filename.variables[var].shape

#show_vars(m1)
#m=m2.variables['CLDICE'][0]-m1.variables['CLDICE'][0]
#print m.max()
#print m1.variables['T'][:]
# print m1.variables['CLOUD'][0][:,40,67]
#print '---------'
#print m1.variables['time']
#for i in range(14,21):
#    x=m1.variables['Z3'][0][i,:,80]
#    plt.plot(lat, x)
#plt.show()
#print Deltah.mean(axis=(1,2))/Deltah.std(axis=(1,2))

# Libradtran input file
def constants(dataset):
    '''
    Returns arrays of latitude, longitude, level and constant P0, which are 
    common throughout datasets.
    '''
    latt, lonn = dataset.variables['lat'][:], dataset.variables['lon'][:]-180.0
    lat, lon = [], []                       # deg, deg
    for lt in latt:
        if lt>=0: lat.append('N '+str(lt))
        elif lt<0: lat.append('S '+str(-lt))
    for ln in lonn:
        if ln>=0: lon.append('E '+str(ln))
        elif ln<0: lon.append('W '+str(-ln))
    lev = dataset.variables['lev'][:]       # level
    P0 = dataset.variables['P0'][:]         # Pa
    # !!!(test time evolution) add 1410 years to fit year between 1900-2000
    date = str(dataset.variables['date'][0]+14100000)
    time = date[:-4].zfill(4), date[-4:-2], date[-2:]
    return P0, lat, lon, lev, time

def variables(dataset, nlat, nlon):
    '''
    Returns standard variables (altitude, pressure, temperature, air density) 
    from given dataset. 
    '''
    P0 = constants(dataset)[0]                                  # Pa
    PS = dataset.variables['PS'][0]                             # Pa
    A = dataset.variables['hyam'][:]                            # fraction
    B = dataset.variables['hybm'][:]                            # fraction
    P = (A[:,None,None]*P0 + B[:,None,None]*PS[None,:,:])/100.0 # mb
    T = dataset.variables['T'][0]                               # K
    air_rho = P/(kb*T)*1.0e-4                                   # cm^-3
    z = dataset.variables['Z3'][0]/1000.0                       # km
    rQ = dataset.variables['RELHUM'][0]                         # percent
    wc = dataset.variables['CLDLIQ'][0]*air0*T0/P0*P/T*100.0    # g cm^-3
    ic = dataset.variables['CLDICE'][0]*air0*T0/P0*P/T*100.0    # g cm^-3
    cfrac = dataset.variables['CLOUD'][0]                       # fraction
    return z[:,nlat,nlon], P[:,nlat,nlon], T[:,nlat,nlon], \
           air_rho[:,nlat,nlon], rQ[:,nlat,nlon], \
           wc[:,nlat,nlon], ic[:,nlat,nlon], cfrac[:,nlat,nlon]

def spec_file(spec, lat, lon, z, Q, u):
    '''
    Creates mol_files for libradtran input.
    '''
    spec_file = open("libradtran/{0}_file.dat".format(spec), 'w')
    spec_file.write('# {2} profile for (lat, lon) = ({0}, {1})\n'
                    .format(lat, lon, spec))
    spec_file.write('# {0:>8} {1:>11}\n'.format('z (km)', spec+' '+u))
    for i in range(len(z)):
        spec_file.write('{0:>10.3f} {1:>11.5f}\n'.format(z[i], Q[i]))
    spec_file.close()

def cloud_file(cf, wc, ic, lat, lon, z):
    '''
    Creates cloud files for libradtran input.
    '''
    wc_file = open("libradtran/wc_file.dat", 'w')
    ic_file = open("libradtran/ic_file.dat", 'w')
    cf_file = open("libradtran/cloudfrac_file.dat", 'w')
    
    wc_file.write('# Water cloud profile for (lat, lon) = ({0}, {1})\n'
                    .format(lat, lon))
    wc_file.write('# {0:>8} {1:>17} {2:>14}\n'
                  .format('z (km)', 'LWC (g cm^-3)', 'R_eff (um)'))
    ic_file.write('# Ice cloud profile for (lat, lon) = ({0}, {1})\n'
                    .format(lat, lon))
    ic_file.write('# {0:>8} {1:>17} {2:>14}\n'
                  .format('z (km)', 'LWC (g cm^-3)', 'R_eff (um)'))
    cf_file.write('# Cloud fraction profile for (lat, lon) = ({0}, {1})\n'
                    .format(lat, lon))
    cf_file.write('# {0:>8} {1:>12}\n'.format('z (km)', 'CF (0-1)'))
    
    for i in range(len(z)):
        wc_file.write('{0:>10.3f} {1:>17.5E} {2:>14.1f}\n'
                  .format(z[i], wc[i], 10.0))
        ic_file.write('{0:>10.3f} {1:>17.5E} {2:>14.1f}\n'
                  .format(z[i], ic[i], 20.0))
        cf_file.write('{0:>10.3f} {1:>12.4E}\n'.format(z[i], cf[i]))
    wc_file.close()
    ic_file.close()
    cf_file.close()

def atmosphere_profile(dataset, nlat, nlon, species=['H2O'], clouds=True):
    '''
    Creates atmosphere file for libradtran input.
    '''
    lat, lon = constants(dataset)[1][nlat], constants(dataset)[2][nlon]
    z, P, T, air_rho, rQ, wc, ic, cf = variables(dataset, nlat, nlon)
    
    atmosphere_file = open("libradtran/atmosphere_file.dat", 'w')
    atmosphere_file.write('# Atmosheric profile for (lat, lon) = ({0}, {1})\n'
                          .format(lat, lon))    
    atmosphere_file.write('# {0:>8} {1:>12} {2:>10} {3:>14}\n'
                          .format('z (km)', 'P (mb)', 'T (K)', 'air (cm^-3)'))
    for i in range(len(z)):
        atmosphere_file.write('{0:>10.3f} {1:>12.5f} {2:>10.3f} {3:>14.5E}\n'
                              .format(z[i], P[i], T[i], air_rho[i]))
    atmosphere_file.close()
    
    for spec in species: # Add more species !!!
        if spec=='H2O': Q, u = rQ, ' (rh)'
        spec_file(spec, lat, lon, z, Q, u)
    
    if clouds: # Fix eff radius !!!. Maybe separate ic and wc
        cloud_file(cf, wc, ic, lat, lon, z)
        

def libra_input(dataset, nlat, nlon, wvl=[2900, 3000], source='solar', 
                param='lowtran', species=['H2O'], clouds=True, 
                error='verbose'):
    '''
    Creates a libradtran input file along with the other necessary files 
    corresponding to the dataset and the (lat, lon) coordinates. 
    Parameter error corresponds to the error handling during the execution of 
    the libradtran code.
    
    source: solar => edir == 0 throughout SW, thermal => edir != 0 throughout LW
    '''
    inp = open("libradtran/test.inp", 'w')
    
    # !!! why spline does not work for thermal? how to make grid of 1nm?
    inp.write('# Wavelength grid and source\n')
    if source=='solar':
        inp.write('wavelength {0} {1} # nm\n'.format(wvl[0], wvl[1]))
        inp.write('source solar ../data/solar_flux/kurudz_1.0nm.dat per_nm\n')
    elif source=='thermal':
        inp.write('wavelength {0} {1} # nm\n'.format(wvl[0], wvl[1]))
        inp.write('source thermal\n')
        #inp.write('spline {0} {1} 1# nm\n'.format(wvl[0], wvl[1]))
    inp.write('\n')
    
    # currently takes standard US profiles of all species except for H2O & CO2
    inp.write('# Atmospheric properties\n')
    atmosphere_profile(dataset, nlat, nlon, species=species, clouds=clouds)
    inp.write('atmosphere_file atmosphere_file.dat\n')
    inp.write('albedo 0.1\n')
    inp.write('## Species\n')
    if 'H2O' in species:
        inp.write('mol_file H2O H2O_file.dat rh\n')
    inp.write('mol_modify CO2 8.45e21 CM_2\n')
    if clouds:
        inp.write('## Clouds\n')
        inp.write('wc_file 1D wc_file.dat\n')
        inp.write('wc_properties mie interpolate\n')
        inp.write('ic_file 1D ic_file.dat\n')
        inp.write('ic_properties yang interpolate\n')
        inp.write('cloud_fraction_file cloudfrac_file.dat\n')
        inp.write('cloud_overlap maxrand\n')
    inp.write('\n')
    
    inp.write('# Geometry\n')
    lat, lon = constants(dataset)[1][nlat], constants(dataset)[2][nlon]
    YYYY, MM, DD = constants(dataset)[-1]
    inp.write('#latitude {0}\n'.format(lat))
    inp.write('#longitude {0}\n'.format(lon))
    inp.write('#time {0} {1} {2} 00 00 00\n'.format(YYYY, MM, DD))
    inp.write('\n')
    
    inp.write('# Radiative transfer equation\n')
    inp.write('rte_solver disort\n')
    inp.write('number_of_streams 6\n')
    inp.write('mol_abs_param {0}\n'.format(param))
    inp.write('\n')
    
    inp.write('# Output\n')
    inp.write('zout atm_levels\n')
    inp.write('output_user wavelength zout T rho_AIR edir edn eup heat\n')
    inp.write('output_process per_nm\n')
    inp.write('\n')
    
    inp.write('{0}\n'.format(error))
    inp.close()

libra_input(m1, 40, 67, wvl=[3000, 10000], source='thermal', 
            param='lowtran', species=['H2O'], clouds=True, 
            error='verbose')







