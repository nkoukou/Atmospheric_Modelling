'''
This module reads netcdf data files from the CCSM3 experiment and creates 
input for the libradtran code.
'''

from netCDF4 import Dataset
import numpy as np
from numpy import sin, cos, tan
import matplotlib.pylab as plt

# Physical constants and pi
pi = np.pi
kb=1.38e-23      # kg m^2 s^-2 K^-1
T0 = 273.15      # K
air0 = 1.2754e-3 # g cm^-3 (air density at T0 and P0=100kPa)

# Import netCDF datasets
def loadnetcdf(month, year, control=True):
    '''
    Loads netCDF history dataset for given year and month. 
    All years have 365 days.
    '''
    if control:
        fdir = "datasets/ccsm/control/b30.031.cam2.h0.{0}-{1}.nc"\
               .format(str(year).zfill(4), str(month).zfill(2))
    else:
        fdir = "datasets/ccsm/double/b30.032a.cam2.h0.{0}-{1}.nc"\
               .format(str(year).zfill(4), str(month).zfill(2))
    return Dataset(fdir)

def avgnetcdf(month, years=[500,550], control=True):
    '''
    Averages all netcdf datasets between year boundaries (years) for given 
    month, producing a netcdf dataset.
    '''
    if control:
        dataset = Dataset("datasets/ccsm/control/avg{0}"
                          .format(str(month).zfill(2)), 'w')
    else:
        dataset = Dataset("datasets/ccsm/double/avg{0}"
                          .format(str(month).zfill(2)), 'w')
    end = loadnetcdf(year=years[1], month=month, control=control)
    count, yrs = 1, []
    for yr in range(years[0], years[1]):
        try: loadnetcdf(year=yr, month=month, control=control)
        except: continue
        count += 1
        yrs.append(yr)
    
    for dim in end.dimensions.keys():
        dataset.createDimension(dim, end.dimensions[dim].size)
    for var in end.variables.keys():
        varr = end.variables[var]
        new = dataset.createVariable(var, varr.dtype, varr.dimensions)
        if not var.isupper(): # !!! add units with try&except
            new[:] = varr[:]
            continue
        new.units = varr.units
        new[:] = varr[:]/count
        for yr in yrs:
            temp = loadnetcdf(year=yr, month=month, control=control)
            new[:] += temp.variables[var][:]/count
    return dataset

def avg_all(years=[500, 550]):
    '''
    Averages all netcdf datasets between year boundaries (years), by running 
    avg_netcdf for all months and both control and double CO2 runs.
    '''
    for month in range(1,13):
        avgnetcdf(month, years=years, control=True)
        avgnetcdf(month, years=years, control=False)
        print 'DONE', month

def select_avg(month, control=True):
    '''
    Selects averaged netcdf file corresponding to the given month from either 
    the control or the double CO2 runs. 
    '''
    if control:
        fdir = "datasets/ccsm/control/avg{0}".format(str(month).zfill(2))
    else:
        fdir = "datasets/ccsm/double/avg{0}".format(str(month).zfill(2))
    return Dataset(fdir)

# Testing data
def show_vars(filename):
    for var in filename.variables.keys():
        print var, filename.variables[var].shape

n=select_avg(4, False)
#show_vars(n)

# Read netCDF datasets
def grid(dataset):
    '''
    Returns arrays of latitude, longitude, level and constant P0, which are 
    common throughout datasets.
    
    lat and lon are converted into lists (each element is of the form N 80 for 
    latitude of 80 degrees North).
    '''
    lat = np.around(dataset.variables['lat'][:], 3) # deg
    lon = np.around(dataset.variables['lon'][:], 3) # deg
    lev = dataset.variables['lev'][:]               # level
    time = dataset.variables['time'][:]%365         # days since start of year
    return lat, lon, lev, time

def variables(dataset, nlat, nlon):
    '''
    Returns standard variables (altitude, pressure, temperature, air density, 
    relative humidity, cloud properties) from given dataset.
    
    Altitude is taken as the geopotential height.
    
    Pressure P is calculated by the atmospheric hybrid sigma coefficients A and 
    B along with the standard and surface pressures P0 and PS respectively with 
    the formula P = A*P0 + B*PS.
    
    Temperature is the radiative temperature.
    
    The number density of air is calculated as P/(kb*T) to preserve 
    hydrodynamic equilibrium.
    
    Cloud liquid and ice condensates are imported in units of fraction =
    = (kg of vapour water)/(kg of dry air). Using the reference air mass 
    density air0 and that P is proportional to rho*T, the mass density of air 
    rho is calculated at all grid points and the units are converted to vapour 
    mass density = fraction*rho. 
    '''
    P0 = dataset.variables['P0'][:]                             # Pa
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

def sza_calc(dataset, nlat):
    '''
    Calculates the solar zenith angle (at a particular grid point specified by 
    latitude and longitude) averaged over a month as an approximation to the 
    monthly averaged data from the CCSM experiments.
    
    Declination dec is calculated throughout the year by taking the reference 
    maximum declination of the Sun as 23.45 and starting from 0.0 on March 21st 
    (day 80 since start of the year).
    
    Standard spherical geometry is used for the rest of the calculations.
    '''
    lat, time = grid(dataset)[0], grid(dataset)[-1]
    lat = np.radians(lat[nlat])
    if time==59.0: dt = 28
    elif time in [120.0, 181.0, 273.0, 334.0]: dt = 30
    else: dt = 31
    time = np.arange(time+1-dt, time+1)
    
    dec = pi * 23.45/180 * sin(2*pi/365 * (time - 80))
    hour_set = np.arccos(-tan(lat)*tan(dec))
    hour_rise = -hour_set
    
    cos_z = 1/(2*pi) * (sin(lat)*sin(dec)*(hour_set-hour_rise) + 
                        cos(lat)*cos(dec)*(sin(hour_set)-sin(hour_rise)))
    z = np.around(180/pi * np.arccos(np.mean(cos_z)), 3)
    return z

# Create libradtran input files
def spec_file(spec, lat, lon, z, Q, u):
    '''
    Creates mol_files for libradtran input. They indicate species densities in 
    the atmosphere.
    
    Currently CO2 and H2O are specified by the model and other gas species are 
    specified from standard US profile datasets.
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
    
    Water and ice cloud droplets have an effective radius of 10 and 20 microns 
    respectively.
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
    lat, lon = grid(dataset)[0][nlat], grid(dataset)[1][nlon]
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
    
    for spec in species:
        if spec=='H2O': Q, u = rQ, ' (rh)'
        spec_file(spec, lat, lon, z, Q, u)
    
    if clouds:
        cloud_file(cf, wc, ic, lat, lon, z)

def libra_input(dataset, nlat, nlon, wvl=[250, 5000], source='solar', 
                species=['H2O'], clouds=True, control=True, error='verbose'):
    '''
    Creates a libradtran input file along with the other necessary files 
    corresponding to the dataset and the (lat, lon) coordinates.
    
    reptran and lowtran are used for solar (SW radiation => edir != 0) and 
    thermal (LW radiation => edir==0) calculations respectively.
    
    species refer to the gas species that are present in the atmosphere. CO2 is 
    set to a specific value throughout the atmosphere (8.45e21 cm^2 for the 
    control run). Data for other species not included in the list are taken 
    from US standard tables.
    
    Parameter error corresponds to the error handling during the execution of 
    the libradtran code.
    
    !!! try kato, do I just double CO2 for experiment?
    '''
    inp = open("libradtran/test.inp", 'w')
    
    inp.write('# Wavelength grid and source\n')
    if source=='solar':
        inp.write('wavelength {0} {1}    # nm\n'.format(wvl[0], wvl[1]))
        inp.write('source solar ../data/solar_flux/kurudz_1.0nm.dat per_nm\n')
        param='reptran'
    elif source=='thermal':
        inp.write('wavelength {0} {1}    # nm\n'.format(wvl[0], wvl[1]))
        inp.write('source thermal\n')
        #inp.write('spline {0} {1} 1    # nm\n'.format(wvl[0], wvl[1]))
        param='lowtran'
    inp.write('\n')
    
    inp.write('# Atmospheric properties\n')
    atmosphere_profile(dataset, nlat, nlon, species=species, clouds=clouds)
    inp.write('atmosphere_file atmosphere_file.dat\n')
    inp.write('albedo 0.1\n')
    inp.write('## Species\n')
    if 'H2O' in species:
        inp.write('mol_file H2O H2O_file.dat rh\n')
    if control:
        CO2_value = 8.45e21     # cm^2
        control = 'Control run'
    else:
        CO2_value = 2 * 8.45e21 # cm^2
        control = 'Double CO2 run'
    inp.write('mol_modify CO2 {0} CM_2    # {1}\n'.format(CO2_value, control))
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
    lat, lon = grid(dataset)[0][nlat], grid(dataset)[1][nlon]
    sza = sza_calc(dataset, nlat)
    inp.write('#latitude {0}    # (-90, 90)\n'.format(lat))
    inp.write('#longitude {0}    # (0, 360)\n'.format(lon))
    inp.write('sza {0}\n'.format(sza))
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

libra_input(n, 30, 60, wvl=[250, 5000], source='solar', species=['H2O'], 
            clouds=True, control=False, error='verbose')



'''
m = loadnetcdf(4,500)
n=select_avg(4,False)
#show_vars(n)
#diff=m500.variables['PS'][0]+m525.variables['PS'][0]+m550['PS'][0]-3*n['PS'][0]
#minn = diff.min()/m550.variables['PS'][0]; maxx=diff.max()/m550.variables['PS'][0]
#print minn.max(), maxx.max()
#print m1.variables['CLOUD'][0][:,40,67]
#print m550.variables['hybm'][:] - m525.variables['hybm'][:]
#print '---------'
#print n.variables['T'][:].shape
#for i in range(14,21):
#    x=m1.variables['Z3'][0][i,:,80]
#    plt.plot(lat, x)
#plt.show()
#print Deltah.mean(axis=(1,2))/Deltah.std(axis=(1,2))
'''




