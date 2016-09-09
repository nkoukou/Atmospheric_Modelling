'''
This module reads netcdf data files from the CCSM experiment.
'''

from netCDF4 import Dataset
import numpy as np

# Import netCDF files
m1 = "datasets/ccsm/b30.031.cam2.h0.0359-01.nc"
m2 = "datasets/ccsm/b30.031.cam2.h0.0359-02.nc"
m1 = Dataset(m1)
m2 = Dataset(m2)

def show_vars(filename):
    for var in filename.variables.keys():
        print var, filename.variables[var].shape

lon = m1.variables['lon'][:]
lat = m1.variables['lat'][:]
lev = m1.variables['lev'][:]
P0 = m1.variables['P0'][:]
PS = m1.variables['PS'][0]
A = m1.variables['hyam'][:]
B = m1.variables['hybm'][:]
P = A[:,None,None]*P0 + B[:,None,None]*PS[None,:,:]
T = m1.variables['T'][0]
TS = m1.variables['TS'][0]
Tsurf = (m1.variables['TSMN'][0]+m1.variables['TSMX'][0])/2
DeltaT = T - TS[None,:,:]

lapse = -0.0065; R = 8.31432; g = 9.80665; M = 0.0289644; c = -R*lapse/g/M
Deltah = (Tsurf/lapse)*((P/PS)**c-1)

#print P.shape
#print '---------'
#print PS
print Deltah[20,20,:]
#print Deltah.mean(axis=(1,2))/Deltah.std(axis=(1,2))
#show_vars(m1)

'''
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
