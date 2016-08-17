'''
This module reads OSSE datasets from '/net/shamal/disk2/truths/' directory.
'''

from scipy.io import readsav
import pickle as pk
import numpy as np

# Latitude and longitude datasets

griddir = "link_osse_sw/2000/wlscale_lat_lon.sav"

coord_sw = readsav(griddir, python_dict=True)
# earth_grid = {'lat':coord_sw['lat'], 'lon':coord_sw['lon']}
a,b = np.meshgrid(coord_sw['lon'],coord_sw['lat'])
earth_grid = np.dstack((b,a))
np.save("datasets/earth_grid",earth_grid)

# Global mean datasets

gm_sw00dir = "link_osse_sw/global_mean/sw_osse_a2_global_mean_2000.sav"
gm_lw00dir = "link_osse_lw/global_mean/lw_osse_a2_global_mean_2000.sav"
gm_lres_sw00 = "link_osse_sw/global_mean/sw_osse_a2_global_mean_lres_2000.sav"
gm_clr_lw00 = "link_osse_lw/global_mean/lw_osse_a2_global_mean_2000_clrsky.sav"
gm_sw99dir = "link_osse_sw/global_mean/sw_osse_a2_global_mean_2099.sav"
gm_lw99dir = "link_osse_lw/global_mean/lw_osse_a2_global_mean_2099.sav"
gm_lres_sw99 = "link_osse_sw/global_mean/sw_osse_a2_global_mean_lres_2099.sav"
gm_clr_lw99 = "link_osse_lw/global_mean/lw_osse_a2_global_mean_2000_clrsky.sav"

swkey, lwkey = 'global_mean_sw_rad', 'global_mean_lw_rad'
lreskey = 'global_mean_sw_rad_lres'

gm_sw00_dict = readsav(gm_sw00dir, python_dict=True)
gm_lw00_dict = readsav(gm_lw00dir, python_dict=True)
gm_lres_sw00_dict = readsav(gm_lres_sw00, python_dict=True)

wvlen = gm_sw00_dict['wavelength']
np.save("datasets/wvlen",wvlen)
wvnum = gm_lw00_dict['wavenumber']
np.save("datasets/wvnum",wvnum)
wvlen_lres = gm_lres_sw00_dict['wavelength_lres']
np.save("datasets/wvlen_lres",wvlen_lres)

gm_sw00 = gm_sw00_dict[swkey]
gm_lw00 = gm_lw00_dict[lwkey] 
gm_lres_sw00 = gm_lres_sw00_dict[lreskey]
gm_clr_lw00 = readsav(gm_clr_lw00, python_dict=True)[lwkey]
gm_sw99 = readsav(gm_sw99dir, python_dict=True)[swkey]
gm_lw99 = readsav(gm_lw99dir, python_dict=True)[lwkey]
gm_lres_sw99 = readsav(gm_lres_sw99, python_dict=True)[lreskey]
gm_clr_lw99 = readsav(gm_clr_lw99, python_dict=True)[lwkey]

np.save("datasets/gm_sw00",gm_sw00)
np.save("datasets/gm_lw00",gm_lw00)
np.save("datasets/gm_lres_sw00",gm_lres_sw00)
np.save("datasets/gm_clr_lw00",gm_clr_lw00)
np.save("datasets/gm_sw99",gm_sw99)
np.save("datasets/gm_lw99",gm_lw99)
np.save("datasets/gm_lres_sw99",gm_lres_sw99)
np.save("datasets/gm_clr_lw99",gm_clr_lw99)

# 2000 & 2099 SW and LW datasets

core = "/b30.042a.cam2.h0.20"
for m in range(1,13):
    if m<10: mstr = '0'+str(m)
    else: mstr = str(m)
    
    sw00dir = "link_osse_sw/2000%s00-%s.sav" % (core, mstr)
    lw00dir = "link_osse_lw/2000%s00-%s.sav" % (core, mstr)
    swl00dir = "link_osse_sw/2000%s00-%s_lres.sav" % (core, mstr)
    lwc00dir = "link_osse_lw/2000%s00-%s_clrsky.sav" % (core, mstr)
    sw99dir = "link_osse_sw/2099%s99-%s.sav" % (core, mstr)
    lw99dir = "link_osse_lw/2099%s99-%s.sav" % (core, mstr)
    swl99dir = "link_osse_sw/2099%s99-%s_lres.sav" % (core, mstr)
    lwc99dir = "link_osse_lw/2099%s99-%s_clrsky.sav" % (core, mstr)
    
    fid = readsav(sw00dir, python_dict=True)['sw_rad']
    name = 'datasets/sw00%s' % (mstr)
    np.save(name,fid)
    fid = readsav(lw00dir, python_dict=True)['lw_rad']
    name = 'datasets/lw00%s' % (mstr)
    np.save(name,fid)
    fid = readsav(swl00dir, python_dict=True)['sw_rad_lres']
    name = 'datasets/lres_sw00%s' % (mstr)
    np.save(name,fid)
    fid = readsav(lwc00dir, python_dict=True)['lw_rad']
    name = 'datasets/clr_lw00%s' % (mstr)
    np.save(name,fid)
    fid = readsav(sw99dir, python_dict=True)['sw_rad']
    name = 'datasets/sw99%s' % (mstr)
    np.save(name,fid)
    fid = readsav(lw99dir, python_dict=True)['lw_rad']
    name = 'datasets/lw99%s' % (mstr)
    np.save(name,fid)
    fid = readsav(swl99dir, python_dict=True)['sw_rad_lres']
    name = 'datasets/lres_sw99%s' % (mstr)
    np.save(name,fid)
    fid = readsav(lwc99dir, python_dict=True)['lw_rad']
    name = 'datasets/clr_lw99%s' % (mstr)
    np.save(name,fid)



