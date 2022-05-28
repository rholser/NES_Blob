%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Update: 25-May-2022

%Imports WOA18 climatology data from netCDF files and saves desired
%measurements as .mat file for ease of use in other scripts.

%This step must be completed prior to Import_CalcAnomalyBinned_MEOP_V2.m

%% Load Climatology Temperature
files=dir('*.nc');

CLat=ncread(files(1).name,'lat');
CLon=ncread(files(1).name,'lon');
CDepth=ncread(files(1).name,'depth');

%Create 4-D array for temperature data so it can be indexed by
%(long,lat,depth,month)
CTemp=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
CTemp2=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
CTsd=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
for i=1:12
    CTemp(:,:,:,i)=ncread(files(i).name,'t_an');
    CTemp2(:,:,:,i)=ncread(files(i).name,'t_mn');
    CTsd(:,:,:,i)=ncread(files(i).name,'t_sd');
end

%% Load Climatology Salinity
files=dir('*.nc');

%Create 4-D array for temperature data so it can be indexed by
%(long,lat,depth,month)
CSal=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
CSal2=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
CSsd=NaN(size(CLon,1),size(CLat,1),size(CDepth,1),12);
for i=1:12
    CSal(:,:,:,i)=ncread(files(i).name,'s_an');
    CSal2(:,:,:,i)=ncread(files(i).name,'s_mn');
    CSsd(:,:,:,i)=ncread(files(i).name,'s_sd');
end
save('D:\Dropbox\MATLAB\Chapter 1\Data\Climatology\Climatology.mat','-v7.3',...
    'CDepth','CLat','CLon','CTemp','CTemp2','CTsd','CSal','CSal2','CSsd');