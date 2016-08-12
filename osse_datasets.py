'''
This module contains datasets which are used in the IDL script 'fentropy.pro'.
'''

from scipy.io import readsav

# OSSE datasets

sw_gm00 = "/net/shamal/disk2/truths/osse_sw/idl_format/global_mean/\
sw_osse_a2_global_mean_2000.sav"

sw_gm99 = "/net/shamal/disk2/truths/osse_sw/idl_format/global_mean/\
sw_osse_a2_global_mean_2099.sav"

lw_gm00 = "/net/shamal/disk2/truths/osse_lw/idl_format/global_mean/\
lw_osse_a2_global_mean_2000.sav"

lw_gm99 = "/net/shamal/disk2/truths/osse_lw/idl_format/global_mean/\
lw_osse_a2_global_mean_2099.sav"


sw_gm00_dict = readsav(sw_gm00, python_dict=True)
sw_gm99_dict = readsav(sw_gm99, python_dict=True)
lw_gm00_dict = readsav(lw_gm00, python_dict=True)
lw_gm99_dict = readsav(lw_gm99, python_dict=True)

'''
# In case there is no gedit display use pickle.
import pickle as pk
nodis=open('nodis.pkl','wb')
pk.dump(sw_gm00_dict,nodis)
pk.dump(sw_gm99_dict,nodis)
pk.dump(lw_gm00_dict,nodis)
pk.dump(lw_gm99_dict,nodis)
nodis.close()
'''



