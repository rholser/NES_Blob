%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 28-May-2022

%%%%%Required scripts, functions, data, and toolboxes:
    %lldistkm 
    %All_CTD_Data_Bin_V2.mat (from Import_CalcAnomalyBinned_MEOP_V2.m)
    %Argo profiling float data from the World Ocean Database (WOD_PFL2) in
        %netCDF format
%%%%%%%

%This script will import seal CTD data and find Argo data that matches given 
%a user-defined time and location window.  
%Use linear models to evaluate the variance between the two data sets and
%plot the results.

%This script will also find Argo casts that are within the same
%temporal-spatiall window and conduct the same comparison within the Argo
%dataset and generate comparable plots

%% Find Matches between MEOP CTD and Argo given time/location window

%This step can take multiple hours

clear
load('All_CTD_Data_V2.mat');
DateConversion=datenum(1770,01,01,0,00,0);
CTD_Data.ArgoDate=CTD_Data.JulDate-DateConversion;

%only need to search for one date/time/location per cast
Subset=CTD_Data(CTD_Data.Depth==1,:); 

%identify Argo data files present
files=dir('*.nc');

%Define window size
days=1;
degs=0.25;

%Preallocate array to receive location matches (rows are seal casts, colums are
%Argo casts).  If match is found, SealxArgo cell will be populated with
%Argo case identifier that can be used to generate the file name containing
%the CTD data for that cast.
Argo_ind=NaN(size(Subset,1),size(files,1));

%Read in Argo index files to search through and load identifiers of Argo
%casts when matches found.
Lat=ncread(files(1).name,'lat');
Lat=[Lat;ncread(files(2).name,'lat')];
Lon=ncread(files(1).name,'lon');
Lon=[Lon;ncread(files(2).name,'lon')];
Date=ncread(files(1).name,'time');
Date=[Date;ncread(files(2).name,'time')];
Cast=ncread(files(1).name,'cast');
Cast=[Cast;ncread(files(2).name,'cast')];

%Search for matches.
tic
for j=1:size(Subset,1)
    for i=1:size(Lat,1)
        if Subset.Lat(j)<=Lat(i)+degs 
            if Subset.Lat(j)>=Lat(i)-degs 
                if Subset.Long(j)<=Lon(i)+degs 
                    if Subset.Long(j)>=Lon(i)-degs 
                        if Subset.ArgoDate(j)<=Date(i)+days
                            if Subset.ArgoDate(j)>=Date(i)-days
                                Argo_ind(j,i)=Cast(i);
                            end
                        end
                    end
                end
            end
        end

    end
end
toc

%Save the resulting array of indices for future work
save('Argo_ind_1day_025deg.mat','-v7.3','Argo_ind')

%% Pair MEOP and Argo T and S Data

%If needed, load MEOP data and matching indices
load('All_CTD_Data_V2.mat');
load('Argo_ind_1day_025deg.mat');

[Seal_ind,Float_ind]=find(Argo_ind>0);

%Identify unique seal and float profiles
Seal_unique=unique(Seal_ind);
Float_unique=unique(Float_ind);

%Create array to populate with data for comparisons
Compare=[];

%For each unique float profile, find the single nearest seal cast that
%matches.  Load the data from each (Temperature and salinity) together into
%a table for analysis.
for i=1:1:size(Float_unique,1)
    
    %load float profile
    ind=find(Float_ind(:,1)==Float_unique(i,1));
    castname=strcat('wod_0',num2str(Argo_ind(Seal_ind(ind(1)),Float_ind(ind(1)))),'O.nc');
    Argo_temp=[ncread(castname,'Temperature') ncread(castname,'Salinity') ncread(castname,'Pressure') ];
    Argo_temp=array2table(Argo_temp,'VariableNames',{'Temp','Sal','Depth'});
    Argo_temp.Lat(:)=ncread(castname,'lat');
    Argo_temp.Lon(:)=ncread(castname,'lon');
    Argo_temp.JulDate(:)=ncread(castname,'time');
    %round depth values to allow matches with seal data
    Argo_temp.DepthR=round(Argo_temp.Depth); 
    
    %If there is only one seal cast matching, load those data
    if size(ind,1)==1
        temp_ind=find(CTD_Data.JulDate==Subset.JulDate(Seal_ind((ind(1)))) & CTD_Data.PTT==Subset.PTT(Seal_ind(ind(1))));
        Seal_temp=CTD_Data(temp_ind(:),:);
        SealTemp=NaN(size(Argo_temp,1),1);
        SealSal=NaN(size(Argo_temp,1),1);
        for j=1:size(Argo_temp,1)
            if Argo_temp.DepthR(j)>0
                try
                    SealTemp(j)=Seal_temp.Temp(Seal_temp.Depth==Argo_temp.DepthR(j));
                    SealSal(j)=Seal_temp.Sal(Seal_temp.Depth==Argo_temp.DepthR(j));
                end
            else
            end
        end
        Comp=[SealTemp SealSal Argo_temp.Temp(:) Argo_temp.Sal(:)];
    %If there are multiple seal casts, choose the one which is closest in
    %space and load it into comparison array.
    else
        for k=1:size(ind,1)
            temp_ind=find(CTD_Data.JulDate==Subset.JulDate(Seal_ind((ind(k)))) & CTD_Data.PTT==Subset.PTT(Seal_ind(ind(k))));
            Locs_temp(k,1)=CTD_Data.Lat(temp_ind(1));
            Locs_temp(k,2)=CTD_Data.Long(temp_ind(1));
            Locs_temp(k,3)=Argo_temp.Lat(1);
            Locs_temp(k,4)=Argo_temp.Lon(1);
            [Locs_temp(k,5) Locs_temp(k,6)]=lldistkm(Locs_temp(k,1:2),Locs_temp(k,3:4));
        end
        %find closest seal cast
        [d l]=min(Locs_temp(:,5)); 
        %find indices of matching seal cast from full data set
        temp_ind=find(CTD_Data.JulDate==Subset.JulDate(Seal_ind((ind(l)))) & CTD_Data.PTT==Subset.PTT(Seal_ind(ind(l))));
        %load seal cast
        Seal_temp=CTD_Data(temp_ind(:),:);
        %Preallocate arrays for seal cast data
        SealTemp=NaN(size(Argo_temp,1),1);
        SealSal=NaN(size(Argo_temp,1),1);
        %Match seal values to depths of Argo cast
        for j=1:size(Argo_temp,1)
            if Argo_temp.DepthR(j)>0
                try
                    SealTemp(j)=Seal_temp.Temp(Seal_temp.Depth==Argo_temp.DepthR(j));
                    SealSal(j)=Seal_temp.Sal(Seal_temp.Depth==Argo_temp.DepthR(j));
                end
            else
            end
        end
        Comp=[SealTemp SealSal Argo_temp.Temp(:) Argo_temp.Sal(:)];
    end
    %Compile comparison data
    Compare=[Compare;Comp];
    clear Argo_temp SealTemp SealSal ind Comp d l Locs_temp castname Seal_temp temp_ind
end
%Conver comparison array into table and name variables
Compare=array2table(Compare,'VariableNames',{'SealTemp','SealSal','ArgoTemp','ArgoSal'});

%% Statistical comparisons
Sind=find(Compare.ArgoSal>=0);

%Fit linear model, SealData~ArgoData
Tmdl=fitlm(Compare.ArgoTemp,Compare.SealTemp);
Tci = coefCI(Tmdl);

Smdl=fitlm(Compare.ArgoSal(Sind),Compare.SealSal(Sind));
Sci = coefCI(Smdl);

%Calculate mean of residuals (abs) as measure of variance in S and T
TResidMean=mean(abs(Tmdl.Residuals.Raw),'omitnan');
SResidMean=mean(abs(Smdl.Residuals.Raw),'omitnan');

TResidStd=std(Tmdl.Residuals.Raw,'omitnan');
SResidStd=std(Smdl.Residuals.Raw,'omitnan');
%Assemble text of linear model equation for legends
testT=strcat('y =', {' '},string(round(Tmdl.Coefficients.Estimate(2),3)),...
    '*x +',{' '},string(round(Tmdl.Coefficients.Estimate(1),3)));
testS=strcat('y =', {' '},string(round(Smdl.Coefficients.Estimate(2),3)),...
    '*x +',{' '},string(round(Smdl.Coefficients.Estimate(1),3)));

%Plot parameters
n=2;
font1=24;
TMin=2;
TMax=24;
SMin=31;
SMax=35;
Font=14;

fig=figure(1);
hold on
p1=plot(Tmdl);
p1(2).Color='k';
p1(2).LineWidth=n;
p1(3).Color='k';
r=refline(1,0);
r.Color=[0.4 0.4 0.4];
r.LineWidth=n;
r.LineStyle='--';
l=legend([p1(2) r],testT,'1:1 Line');
ax=gca;
ax.FontSize=Font;
ax.XLim=[TMin,TMax];
ax.YLim=[TMin,TMax];
ax.XLabel.String='Argo Temperature';
ax.YLabel.String='Seal Temperature';

fig=figure(2);
hold on
p2=plot(Smdl)
p2(2).Color='k';
p2(2).LineWidth=n;
p2(3).Color='k';
r=refline(1,0);
r.Color=[0.4 0.4 0.4];
r.LineWidth=n;
r.LineStyle='--';
l=legend([p2(2) r],testS,'1:1 Line');
ax=gca;
ax.FontSize=Font;
ax.XLim=[SMin,SMax];
ax.YLim=[SMin,SMax];
ax.XLabel.String='Argo Salinity';
ax.YLabel.String='Seal Salinity';

%% Locations of matching casts and map
Seal_unique=unique(Seal_ind);
Float_unique=unique(Float_ind);

Seal_Locations=NaN(size(Seal_unique,1),3);
Float_Locations=NaN(size(Float_unique,1),3);

for i=1:size(Float_unique)
    ind=find(Float_ind==Float_unique(i));
    castname=strcat('wod_0',num2str(Argo_ind(Seal_ind(ind(1)),Float_ind(ind(1)))),'O.nc');
    info=ncinfo(castname);
    Float_Locations(i,1)=ncread(castname,'lat');
    Float_Locations(i,2)=ncread(castname,'lon');
    Float_Locations(i,3)=ncread(castname,'time');

    clear castname ind
end
ind=find(Float_Locations(:,2)<0);
Float_Locations(ind,2)=Float_Locations(ind,2)+360;

for i=1:size(Seal_unique)
    ind=find(Seal_ind==Seal_unique(i));
    temp_ind=find(CTD_Data.JulDate==Subset.JulDate(Seal_ind(ind(1))) & CTD_Data.PTT==Subset.PTT(Seal_ind(ind(1))));
    Seal_temp=CTD_Data(temp_ind(:),:);  
    Seal_Locations(i,1)=Seal_temp.Lat(1);
    Seal_Locations(i,2)=Seal_temp.Long360(1);
    Seal_Locations(i,3)=Seal_temp.ArgoDate(1);

    clear Seal_temp ind temp_ind
end

fig=figure(3);
fig.Color=[1,1,1];
g1=geoscatter(Seal_Locations(:,1),Seal_Locations(:,2),20,'filled');
hold on
g2=geoscatter(Float_Locations(:,1),Float_Locations(:,2),20,'filled');
geobasemap landcover
ax=gca;
hold off
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.Box='on';
ax.LineWidth=1.5;
l=legend([g1 g2],'Seal CTD Casts','Argo Float Casts');


%% Find Argo-Argo matches given time/location window

%Define window size
days=1;
degs=0.25;

%Preallocate array to receive location matches (rows are casts, colums are
%Argo casts).  If match is found, ArgoxArgo cell will be populated with
%Argo case identifier that can be used to generate the file name containing
%the CTD data for that cast.
Argo2_ind=NaN(size(files,1),size(files,1));

%Read in Argo index files to search through and load identifiers of Argo
%casts when matches found.
Lat=ncread(files(1).name,'lat');
Lat=[Lat;ncread(files(2).name,'lat')];
Lon=ncread(files(1).name,'lon');
Lon=[Lon;ncread(files(2).name,'lon')];
Date=ncread(files(1).name,'time');
Date=[Date;ncread(files(2).name,'time')];
Cast=ncread(files(1).name,'cast');
Cast=[Cast;ncread(files(2).name,'cast')];

%Search for matches.
tic
for j=1:size(Lat,1)
    for i=1:size(Lat,1)
        if i==j
            Argo_ind(j,i)=NaN;
        elseif Lat(j)<=Lat(i)+degs 
            if Lat(j)>=Lat(i)-degs 
                if Lon(j)<=Lon(i)+degs 
                    if Lon(j)>=Lon(i)-degs 
                        if Date(j)<=Date(i)+days
                            if Date(j)>=Date(i)-days
                                Argo2_ind(j,i)=Cast(i);
                            end
                        end
                    end
                end
            end
        end

    end
end
toc

%% Argo-Argo Comparison

[Float1_ind,Float2_ind]=find(Argo2_ind>0);

%Unique seal and float profiles
Float1_unique=unique(Float1_ind);
Float2_unique=unique(Float2_ind);

%Create array to populate with data for comparisons
Compare=[];

%For each unique float profile, find the single nearest seal cast that
%matches.  Load the data from each (Temperature and salinity) together into
%a table for analysis.
for i=1:4:344
    
    %load float profile
    castname1=strcat('wod_0',num2str(Argo2_ind(Float1_ind(i),Float2_ind(i))),'O.nc');
    Argo_temp=[ncread(castname1,'Temperature') ncread(castname1,'Salinity') ncread(castname1,'Pressure') ];
    Argo_temp=array2table(Argo_temp,'VariableNames',{'Temp','Sal','Depth'});
    Argo_temp.Lat(:)=ncread(castname1,'lat');
    Argo_temp.Lon(:)=ncread(castname1,'lon');
    Argo_temp.JulDate(:)=ncread(castname1,'time');
    %round depth values to allow matches with seal data
    Argo_temp.DepthR=round(Argo_temp.Depth); 
    
    %load float profile 2
    castname2=strcat('wod_0',num2str(Argo2_ind(Float1_ind(i+1),Float2_ind(i+1))),'O.nc');
    Argo_temp2=[ncread(castname2,'Temperature') ncread(castname2,'Salinity') ncread(castname2,'Pressure') ];
    Argo_temp2=array2table(Argo_temp2,'VariableNames',{'Temp','Sal','Depth'});
    Argo_temp2.Lat(:)=ncread(castname2,'lat');
    Argo_temp2.Lon(:)=ncread(castname2,'lon');
    Argo_temp2.JulDate(:)=ncread(castname2,'time');
    %round depth values to allow matches with seal data
    Argo_temp2.DepthR=round(Argo_temp2.Depth); 
    
    %If there is only one seal cast matching, load those data
    Temp=NaN(size(Argo_temp,1),1);
    Sal=NaN(size(Argo_temp,1),1);
    for j=1:size(Argo_temp,1)
        if Argo_temp.DepthR(j)>0
            try
                Temp(j)=Argo_temp2.Temp(Argo_temp2.DepthR==Argo_temp.DepthR(j));
                Sal(j)=Argo_temp2.Sal(Argo_temp2.DepthR==Argo_temp.DepthR(j));
            end
        else
        end
    end

    Comp=[Temp Sal Argo_temp.Temp(:) Argo_temp.Sal(:)];
    %Compile comparison data
    Compare=[Compare;Comp];
    clear Argo_temp SealTemp SealSal ind Comp d l Locs_temp castname Seal_temp temp_ind
end
%Conver comparison array into table and name variables
Compare=array2table(Compare,'VariableNames',{'Argo2Temp','Argo2Sal','ArgoTemp','ArgoSal'});

%% Argo-Argo Statistical comparisons
Compare(isnan(Compare.Argo2Temp),:)=[];
Sind=find(Compare.ArgoSal<1);
Compare(Sind,:)=[];

%Fit linear model, Argo2Data~ArgoData
Tmdl=fitlm(Compare.ArgoTemp,Compare.Argo2Temp);
Tci = coefCI(Tmdl);

Smdl=fitlm(Compare.ArgoSal,Compare.Argo2Sal);
Sci = coefCI(Smdl);

%Calculate mean of residuals (abs) as measure of variance in S and T
TResidMean=mean(abs(Tmdl.Residuals.Raw),'omitnan');
SResidMean=mean(abs(Smdl.Residuals.Raw),'omitnan');

TResidStd=std(Tmdl.Residuals.Raw,'omitnan');
SResidStd=std(Smdl.Residuals.Raw,'omitnan');
%Assemble text of linear model equation for legends
testT=strcat('y =', {' '},string(round(Tmdl.Coefficients.Estimate(2),3)),...
    '*x +',{' '},string(round(Tmdl.Coefficients.Estimate(1),3)));
testS=strcat('y =', {' '},string(round(Smdl.Coefficients.Estimate(2),3)),...
    '*x +',{' '},string(round(Smdl.Coefficients.Estimate(1),3)));

%Plot parameters
n=2;
font1=24;
TMin=2;
TMax=24;
SMin=31;
SMax=35;
Font=14;

fig=figure(1);
hold on
p1=plot(Tmdl);
p1(2).Color='k';
p1(2).LineWidth=n;
p1(3).Color='k';
r=refline(1,0);
r.Color=[0.4 0.4 0.4];
r.LineWidth=n;
r.LineStyle='--';
l=legend([p1(2) r],testT,'1:1 Line');
ax=gca;
ax.FontSize=Font;
ax.XLim=[TMin,TMax];
ax.YLim=[TMin,TMax];
ax.XLabel.String='Argo Temperature';
ax.YLabel.String='Argo2 Temperature';

fig=figure(2);
hold on
p2=plot(Smdl)
p2(2).Color='k';
p2(2).LineWidth=n;
p2(3).Color='k';
r=refline(1,0);
r.Color=[0.4 0.4 0.4];
r.LineWidth=n;
r.LineStyle='--';
l=legend([p2(2) r],testS,'1:1 Line');
ax=gca;
ax.FontSize=Font;
ax.XLim=[SMin,SMax];
ax.YLim=[SMin,SMax];
ax.XLabel.String='Argo Salinity';
ax.YLabel.String='Argo2 Salinity';