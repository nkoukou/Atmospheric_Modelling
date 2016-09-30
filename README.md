# Entropy and radiation intensity calculations on OSSE and CCSM data

This project was undertaken as a supervised summer internship and is divided 
into two parts. The general purpose of the project is the analysis of radiation 
intensity data obtained by an Observational System Simulation Experiment (OSSE).
The OSSE studied is the Community Climate System Model version 3 (CCSM3). The 
data used from the model are monthly averages and can be accessed at [the 
website of the National Center for Atmospheric Research (NCAR)]
(https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm.output.html).

The first part focused on CCSM3 run b30.042a. The radiative and entropic net 
fluxes that shape Earth's entropy budget were analysed as a function of 
latitude and longitude. Global maps illustrating the net flux were produced, 
while a separate analysis of global mean data was also performed. The relevant 
python modules developed are 1-4 as listed below.

The second part focused on CCSM3 runs b30.031 and b30.032a which are a control 
run and a run with double CO_2 concentration respectively. An atmospheric model 
was produced at a 3D grid (latitude, longitude, altitude) based on data from 
the runs. The model considers solar radiation, pressure, temperature, water 
vapour, cloud properties etc. and was used as the input of the [libradtran 
radiation code](www.libradtran.org/doku.php) which calculates direct and 
diffused fluxes based on the atmosheric profile. Finally, analysis was 
performed based on the libradtran output fluxes and, in particular, estimates 
of the total material entropy of the Earth were calculated. The calculations 
can be used to investigate the impact of CO_2 on Earth's entropy. The relevant 
python modules developed are 5-6 as listed below.

-------
## Modules

Descriptions can be found at the top of each module.

1. osse_datasets.py
2. ent_datasets.py
3. entropy.py
4. glob_entropy.py
5. ccsm_datasets.py
6. mat_ent.py

-------
## Important notes (listed by module)

1. Data for years 2000 and 2099 are imported.
   
   This module exports data.

2. A factor of 1000 corrects SW (non lres) radiation data.
   Entropy conversion is more accurate with Taylor approximation 
   for low radiation.
   Errors messages have been supressed for entropy calculations because the 
   numpy.where conditionals handle them successfully.
   
   This module exports data.

3. Adapts entropy.pro into python with functions grouping the various
   plots.
   Comments on entropy and flux from above module apply here too (no lres data 
   were used here).
   
   The plotting functions are the main functions of the module.

4. Analyses global data of radiation exported from ent_datasets.py.
   Map colorbars do not have an absolute colour scale throughout all maps.
   
   analyse\_month(), spec\_analysis() and analyse\_month_mini() are the main 
   functions of this module. They print global maps and can also return 
   information on the data used to plot the maps.

5. Creates atmospheric profiles as input for libradtran code. Gets data as 
   netcdf files from the CCSM control and CO_2 doubling experiments.
   Methods of calculation for various atmospheric parameters are explained in 
   the code.
   Solar zenith angle function is unreliable for extreme latitudes (outside the 
   range of (-80, 80) degrees).
   
   libra_input() is the main function of the module. It creates all necessary 
   .inp files.

6. Calculates global material entropy. Gets data from libradtran code output.
   Care should be taken with the plotting functions when used outside the 
   namespace of the flux_output() function. Redesign of the plot functions is 
   encouraged if considered useful on their own.
   
   flux_output() is the main function of the module. It creates two .txt files 
   with calculations and several plots.

Examples of important functions are commented out at the bottom of each script. 
Function arguments are explained in the documentation of the function.

-------
ToDo:

1. Tidy up directories, fix ALL links in code
2. Module 4: Re-export graphs (colorbar consistency, framing)
3. Module 5: Check time evolution between years 500 and 550, add more years for 
   averaging
4. Modules 5,6: write bash script to run code for many lat & lon & months and 
   run the libradtran code in between. This would help automate material entropy
   calculations.
