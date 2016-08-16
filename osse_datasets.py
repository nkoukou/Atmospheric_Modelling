'''
This module reads OSSE datasets from '/net/shamal/disk2/truths/' directory.
'''

from scipy.io import readsav
import numpy as np

# Latitude and longitude datasets

griddir = "link_osse_sw/2000/wlscale_lat_lon.sav"

coord_sw = readsav(griddir, python_dict=True)
earth_grid = {'lat':coord_sw['lat'], 'lon':coord_sw['lon']}

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

years = ['00','99']
core = "/b30.042a.cam2.h0.20"
sw_dict = {}
lw_dict = {}
lres_sw_dict = {}
clr_lw_dict = {}
for ystr in years:
    for m in range(1,13):
        if m<10: mstr = '0'+str(m)
        else: mstr = str(m)
        swdir = "link_osse_sw/20"+ystr+core+ystr+"-"+mstr+".sav"
        lwdir = "link_osse_lw/20"+ystr+core+ystr+"-"+mstr+".sav"
        swldir = "link_osse_sw/20"+ystr+core+ystr+"-"+mstr+"_lres.sav"
        lwcdir = "link_osse_lw/20"+ystr+core+ystr+"-"+mstr+"_clrsky.sav"
        sw_dict['sw'+ystr+mstr] = readsav(swdir, python_dict=True)['sw_rad']
        lw_dict['lw'+ystr+mstr] = readsav(lwdir, python_dict=True)['lw_rad']
        lres_sw_dict['lres_sw'+ystr+mstr] = readsav(swldir,
          python_dict=True)['sw_rad_lres']
        clr_lw_dict['clr_lw'+ystr+mstr] = readsav(lwcdir,
          python_dict=True)['lw_rad']

'''
# 2000 & 2099 clrsky & lres datasets

lres_sw_dict = {}
clr_lw_dict = {}
for ystr in years:
    for m in range(1,13):
        if m<10: mstr = '0'+str(m)
        else: mstr = str(m)
        swldir = "link_osse_sw/20"+ystr+core+ystr+"-"+mstr+"_lres.sav"
        lwcdir = "link_osse_lw/20"+ystr+core+ystr+"-"+mstr+"_clrsky.sav"
        lres_sw_dict['lres_sw'+ystr+mstr+'_dict'] = readsav(swldir,
          python_dict=True)['sw_rad_lres']
        clr_lw_dict['clr_lw'+ystr+mstr+'_dict'] = readsav(lwcdir,
          python_dict=True)['lw_rad']
'''

gm_dict = {'sw00':gm_sw00, 'lw00':gm_lw00, 'sw99':gm_sw99, 'lw99':gm_lw99,
           'wvl':wvlen, 'wvn':wvnum}
lc_gm_dict = {'sw00':gm_lres_sw00, 'lw00':gm_clr_lw00, 'sw99':gm_lres_sw99,
              'lw99':gm_clr_lw99, 'wvl':wvlen_lres, 'wvn':wvnum}






