'''
This module reads netcdf data files from the CCSM3 experiment and creates 
input for the libradtran code.
'''

from netCDF4 import Dataset
import numpy as np
import matplotlib.pylab as plt

kb=1.38e-23 # kg m^2 s^-2 K^-1
# Import netCDF files
m1 = "datasets/ccsm/b30.031.cam2.h0.0359-01.nc"
m2 = "datasets/ccsm/b30.031.cam2.h0.0359-02.nc"
m1 = Dataset(m1)
m2 = Dataset(m2)

def show_vars(filename):
    for var in filename.variables.keys():
        print var, filename.variables[var].shape

lon = m1.variables['lon'][:]                                # deg
lat = m1.variables['lat'][:]                                # deg
lev = m1.variables['lev'][:]                                # level
P0 = m1.variables['P0'][:]                                  # Pa
PS = m1.variables['PS'][0]                                  # Pa
A = m1.variables['hyam'][:]                                 # fraction
B = m1.variables['hybm'][:]                                 # fraction
P = A[:,None,None]*P0 + B[:,None,None]*PS[None,:,:]         # Pa
T = m1.variables['T'][0]                                    # K
TS = m1.variables['TS'][0]                                  # K
air_rho = P/(kb*T)*1.0e-6                                   # cm^-3

#show_vars(m1)
print m1.variables['time']
#print '---------'
#print m1.variables['time']
#for i in range(14,21):
#    x=m1.variables['Z3'][0][i,:,80]
#    plt.plot(lat, x)
#plt.show()
#print Deltah.mean(axis=(1,2))/Deltah.std(axis=(1,2))

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
    date = str(dataset.variables['date'][0]+16000000)
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
    return z[:,nlat,nlon], P[:,nlat,nlon], T[:,nlat,nlon], \
           air_rho[:,nlat,nlon], rQ[:,nlat,nlon]

def atmosphere_profile(dataset, nlat, nlon, species=['H2O']):
    '''
    Creates atmosphere file for libradtran input.
    '''
    lat, lon = constants(dataset)[1][nlat], constants(dataset)[2][nlon]
    z, P, T, air_rho, rQ = variables(dataset, nlat, nlon)
    atmosphere_file = open("libradtran/atmosphere_file.dat", 'w')
    atmosphere_file.write('# Atmosheric profile for (lat, lon) = ({0}, {1})\n'
                          .format(lat, lon))    
    atmosphere_file.write('# {0:>8} {1:>12} {2:>10} {3:>14}\n'
                          .format('z (km)', 'P (mb)', 'T (K)', 'air (cm^-3)'))
    for i in range(len(z)):
        atmosphere_file.write('{0:>10.3f} {1:>12.5f} {2:>10.3f} {3:>14.5E}\n'
                              .format(z[i], P[i], T[i], air_rho[i]))
    atmosphere_file.close()
    
    for spec in species: # !!! rh and rQ must be changed in the loop
        spec_file = open("libradtran/{0}_file.dat".format(spec), 'w')
        spec_file.write('# {2} profile for (lat, lon) = ({0}, {1})\n'
                       .format(lat, lon, spec))
        spec_file.write('# {0:>8} {1:>11}\n'.format('z (km)', spec+' (rh)'))
        for i in range(len(z)):
            spec_file.write('{0:>10.3f} {1:>11.5f}\n'.format(z[i], rQ[i]))
        spec_file.close()

def libra_input(dataset, nlat, nlon, species=['H2O'], error='verbose'):
    '''
    Creates a libradtran input file along with the other necessary files 
    corresponding to the dataset and the (lat, lon) coordinates. 
    Parameter error corresponds to the error handling during the execution of 
    the libradtran code.
    '''
    inp = open("libradtran/test.inp", 'w')
    
    # How to determine wavelength? !!!
    # Add clouds!!!
    inp.write('# Wavelength grid and source\n')
    inp.write('wavelength 250 10000 # nm\n')
    inp.write('source thermal\n')
    inp.write('\n')
    
    # mol_modify and mol_file for all species densities !!!
    inp.write('# Atmospheric properties\n')
    atmosphere_profile(dataset, nlat, nlon, species=species)
    inp.write('atmosphere_file atmosphere_file.dat\n')
    inp.write('mol_file H2O H2O_file.dat rh\n')
    inp.write('mol_modify CO2 8.45e21 CM_2\n')
    inp.write('albedo 0.1\n')
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
    inp.write('mol_abs_param reptran fine\n')
    inp.write('\n')
    
    inp.write('# Output\n')
    inp.write('zout atm_levels\n')
    inp.write('output_user wavelength zout T rho_AIR edir edn eup heat\n')
    inp.write('output_process per_nm\n')
    inp.write('\n')
    
    inp.write('{0}\n'.format(error))
    inp.close()

libra_input(m1, 40, 67)







'''
Tsurf = (m1.variables['TSMN'][0]+m1.variables['TSMX'][0])/2 # K
DeltaT = Tsurf - TS[None,:,:]
lapse = -0.0065; R = 8.31432; g = 9.80665; M = 0.0289644; c = -R*lapse/g/M
Deltah = (Tsurf/lapse)*((P/PS)**c-1)

liq = "datasets/ccsm/b30.031.cam2.h0.CLDLIQ_zav.0200-01_cat_0299-12.nc"
ice = "datasets/ccsm/b30.031.cam2.h0.CLDICE_zav.0200-01_cat_0299-12.nc"
ps = "datasets/ccsm/b30.031.cam2.h0.PS.0200-01_cat_0299-12.nc"
hgh = "datasets/ccsm/test/b30.031.cam2.h0.CLDHGH.0200-01_cat_0299-12.nc"
low = "datasets/ccsm/test/b30.032a.cam2.h0.CLDLOW.0479-01_cat_0578-12.nc"
frac = "datasets/ccsm/b30.031.cam2.h0.CLOUD_zav.0200-01_cat_0299-12.nc"
hum = "datasets/ccsm/b30.031.cam2.h0.Q_zav.0200-01_cat_0299-12.nc"
temp = "datasets/ccsm/b30.031.cam2.h0.T_zav.0200-01_cat_0299-12.nc"
hum2 = "datasets/ccsm/test/b30.032a.cam2.h0.Q_zav.0479-01_cat_0578-12.nc"
temp2 = "datasets/ccsm/test/b30.032a.cam2.h0.T_zav.0479-01_cat_0578-12.nc"
tsoi = "datasets/ccsm/test/b30.040e.clm2.h0.TSOI.2000-01_cat_2099-12.nc"
temp = "datasets/ccsm/test/1999/b30.030e.cam2.h3.SA_Q.1999-01-01_cat_1999-12-31.nc"
temp2 = "datasets/ccsm/test/1999/b30.030e.cam2.h3.Q.1999-11-01_cat_1999-11-30.nc"

liq = Dataset(liq)
ice = Dataset(ice)
frac = Dataset(frac)
hgh = Dataset(hgh)
low = Dataset(low)
hum = Dataset(hum2)
temp = Dataset(temp2)
ps = Dataset(ps)
tsoi = Dataset(tsoi)
temp = Dataset(temp)
temp2 = Dataset(temp2)

#print ice.variables['P0'][:]
#print ice.variables['lev'][:]
#print ice.variables['CLDICE']
#print ice.variables['hyai']
#print ice.variables['hyam'][:]
#print ice.variables['hybi']
#print ice.variables['hybm'][:]
#print ice.variables['ilev']
'''
