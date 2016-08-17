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
wvnum = gm_lw00_dict['wavenumber']
wvlen_lres = gm_lres_sw00_dict['wavelength_lres']

gm_sw00 = gm_sw00_dict[swkey]
gm_lw00 = gm_lw00_dict[lwkey] 
gm_lres_sw00 = gm_lres_sw00_dict[lreskey]
gm_clr_lw00 = readsav(gm_clr_lw00, python_dict=True)[lwkey]
gm_sw99 = readsav(gm_sw99dir, python_dict=True)[swkey]
gm_lw99 = readsav(gm_lw99dir, python_dict=True)[lwkey]
gm_lres_sw99 = readsav(gm_lres_sw99, python_dict=True)[lreskey]
gm_clr_lw99 = readsav(gm_clr_lw99, python_dict=True)[lwkey]

# 2000 & 2099 SW and LW datasets

core = "/b30.042a.cam2.h0.20"
sw00_dict = {}; sw99_dict = {}
lw00_dict = {}; lw99_dict = {}
lres_sw00_dict = {}; lres_sw99_dict = {}
clr_lw00_dict = {}; clr_lw99_dict = {}
for m in range(1,13):
    if m<10: mstr = '0'+str(m)
    else: mstr = str(m)
    sw00dir = "link_osse_sw/2000"+core+"00-"+mstr+".sav"
    lw00dir = "link_osse_lw/2000"+core+"00-"+mstr+".sav"
    swl00dir = "link_osse_sw/2000"+core+"00-"+mstr+"_lres.sav"
    lwc00dir = "link_osse_lw/2000"+core+"00-"+mstr+"_clrsky.sav"
    sw99dir = "link_osse_sw/2099"+core+"99-"+mstr+".sav"
    lw99dir = "link_osse_lw/2099"+core+"99-"+mstr+".sav"
    swl99dir = "link_osse_sw/2099"+core+"99-"+mstr+"_lres.sav"
    lwc99dir = "link_osse_lw/2099"+core+"99-"+mstr+"_clrsky.sav"
    sw00_dict[mstr] = readsav(sw00dir, python_dict=True)['sw_rad']
    lw00_dict[mstr] = readsav(lw00dir, python_dict=True)['lw_rad']
    lres_sw00_dict[mstr] = readsav(swl00dir, python_dict=True)['sw_rad_lres']
    clr_lw00_dict[mstr] = readsav(lwc00dir, python_dict=True)['lw_rad']
    sw99_dict[mstr] = readsav(sw99dir, python_dict=True)['sw_rad']
    lw99_dict[mstr] = readsav(lw99dir, python_dict=True)['lw_rad']
    lres_sw99_dict[mstr] = readsav(swl99dir, python_dict=True)['sw_rad_lres']
    clr_lw99_dict[mstr] = readsav(lwc99dir, python_dict=True)['lw_rad']


gm_dict = {'sw00':gm_sw00, 'lw00':gm_lw00, 'sw99':gm_sw99, 'lw99':gm_lw99,
           'wvl':wvlen, 'wvn':wvnum}
lc_gm_dict = {'sw00':gm_lres_sw00, 'lw00':gm_clr_lw00, 'sw99':gm_lres_sw99,
              'lw99':gm_clr_lw99, 'wvl':wvlen_lres, 'wvn':wvnum}
w00_dict = {'sw':sw00_dict, 'lw':lw00_dict, 'wvl':wvlen, 'wvn':wvnum,
            'grid':earth_grid}
w99_dict = {'sw':sw99_dict, 'lw':lw99_dict, 'wvl':wvlen, 'wvn':wvnum,
            'grid':earth_grid}
lc_w00_dict = {'sw':lres_sw00_dict, 'lw':clr_lw00_dict, 'wvl':wvlen_lres,
               'wvn':wvnum, 'grid':earth_grid}
lc_w99_dict = {'sw':lres_sw99_dict, 'lw':clr_lw99_dict, 'wvl':wvlen_lres,
               'wvn':wvnum, 'grid':earth_grid}

gm_pk = open('gm_dict.pkl','wb')
pk.dump(gm_dict, gm_pk)
gm_pk.close()
lc_gm_pk = open('lc_gm_dict.pkl','wb')
pk.dump(lc_gm_dict, lc_gm_pk)
lc_gm_pk.close()
w00_pk = open('w00_dict.pkl','wb')
pk.dump(w00_dict, w00_pk)
w00_pk.close()
w99_pk = open('w99_dict.pkl','wb')
pk.dump(w99_dict, w99_pk)
w99_pk.close()
lc_w00_pk = open('lc_w00_dict.pkl','wb')
pk.dump(lc_w00_dict, lc_w00_pk)
lc_w00_pk.close()
lc_w99_pk = open('lc_w99_dict.pkl','wb')
pk.dump(lc_w99_dict, lc_w99_pk)
lc_w99_pk.close()







