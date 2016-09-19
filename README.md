Entropy and intensity calculations

Temporarily, issues are listed below:
-------------------------------------
papers / GIT (priority) / bash scripting
-------------------------------------
ToDo:

-1. Request one overlapping years
-2. Lowtran works with thermal, try for solar (vs reptran)
-3. Keep notes and comments of all assumptions and calculations
4. Calculate sza theoretically
5. mat_ent.py (deriv function, calculations, plots, output files)
6. Request many overlapping years (be careful with size)
   Check time evolution
   Make the necessary 3D fields (make netcdf file out of netcdf avging)
7. exported graphs (colorbar consistency, framing)
-------------------------------------
Modules:

1. osse_datasets.py
2. ent_datasets.py
3. entropy.py
4. glob_entropy.py
5. ccsm_datasets.py
6. mat_ent.py
-------------------------------------
Notes (listed by module):

1. Years 2000 and 2099 are imported.

2. A factor of 2000 corrects SW (non lres) radiation data.
   Entropy conversion is more accurate with Taylor approximation 
   for low radiation.
   Errors messages have been supressed for entropy calculations because the 
   numpy.where conditionals handle them successfully.
   Calculated values preserve the accuracy of the physical constants quoted.
   Numerical integration can become more efficient by use of numpy.trapz

3. Adapts entropy.pro into python with functions grouping the various
   plots.
   Comments on entropy and flux carry on from 2. (no lres data were used here)

4. Analyses global data of radiation exported from ent_datasets.py.
   Map colorbars do not have an absolute colour scale which would be the same 
   throughout all maps. 

5. Creates atmospheric profiles as input for libradtran code. Gets data as 
   netcdf files from the CCSM control and CO_2 doubling experiments.
   Methods of calculation for various atmospheric parameters are explained in 
   the code.

6. Calculates global material entropy. Gets data from libradtran code output.
