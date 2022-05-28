# NES_Blob
Code used to import, process, analyze, and plot northern elephant seal-collected oceanographic data compared to World Ocean Atlas 2018 climatologies during the 2014/2015 marine heatwave.

All raw data are publicly available:

Northern elephant seal-collected CTD data can be accessed through the MEOP data portal: https://www.meop.net/database/meop-databases/  or through SEANOE:  https://www.seanoe.org/data/00343/45461/.  Data are stored in netCDF format.

Climatological data used were from World Ocean Atlas 2018 (1981-2010 climatology on a 1x1 degree grid): https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/  Code is written to deal with netCDF format.

All scripts are written in Matlab.

Requires Gibbs Seawater Toolbox (TEOS 10): https://www.teos-10.org/software.htm
Requires sw_dens function: https://github.com/ashao/matlab/blob/master/external/seawater/sw_dens.m

Step 1: run Import_Climatology_Data.m <-- pulls climatology data into .mat file for easy access in later steps
Step 2: run Import_CalcAnomalyBinned_MEOP_V2.m  <-- loads seal CTD data, completed additional quality control and filtering, matches data to climatology from Step 1, and calculates temperature, salinity and spice anomalies.  Saves resulting array as .mat file for later use.

Once data are imported, remaining scripts can be completed in any order:
    Create various plots of anomaly data using 'plot_' scripts
    Calculate water movement along isopycnals between climatological and marine heatwave conditions using an objective function  (Analysis_DensitySurface_Distance_V7.m)
    Calculate mixed layer depth
    Co-locate Argo profiling float data with NES CTD data in space and time; compare measurements
