%Rachel Holser (rholser@ucsc.edu)
%Last Updated: 28-May-2022

%%%%% Required scripts, functions, and toolboxes:
    %GSW toolbox
    %Climatology.mat (from Import_Climatology_Data.m)
%%%%%%%

%Step 1: Import CTD Data
    %Import MEOP QC'd CTD data from netCDF (only retains quality 1, 2, and 3 
    %measurements for temperature and salinity)

    %Calculates Conservative Temperature (using GSW Toolbox)

    %Calculates mean temperatures in depth bins to match Climatology depth resolution.

%Step 2: Load Climatology Data
    %Imported and compiled into .mat file in Import_Climatology_Data.m

%Step 3: Match Climatology Data
    %Use knnsearch to find nearest match in lat, lon, and depth
    %Calculates T and S anomaly from both an and mn for comparison.

%Step 4: Filter Outlier Data 
    %Filters data that are +/- 8 s.d. from mean of complete data set at three 
    %depth segments(<=100m, 100-500m, >=500m)

%Step 5: Calculate sigmaTa for each measurement

%Step 6: Calculate spiciness and spice anomaly

%Step 7: Save data as .mat file


%% Step 1: Import CTD Data
clear
files=dir('*.nc');
CTD_Data_Bin=[];

%load each data file and extract T, S, depth, Lat, Long, JulDate.  Compile
%into one n x 6 table for further analysis
for i=1:size(files,1)
    %rows are depth, columns are casts
    info=ncinfo(files(i).name);
    ID=string(info.Attributes(45).Value);
    ID2=str2num(string(info.Attributes(46).Value));
    Temp=ncread(files(i).name,'TEMP_ADJUSTED');
    Temp_QC=ncread(files(i).name,'TEMP_ADJUSTED_QC');
    Sal=ncread(files(i).name,'PSAL_ADJUSTED');
    Sal_QC=ncread(files(i).name,'PSAL_ADJUSTED_QC');
    Pres=ncread(files(i).name,'PRES_ADJUSTED');
    Pres_QC=ncread(files(i).name,'PRES_ADJUSTED_QC');
    
    %each row is value for cast
    Lat=ncread(files(i).name,'LATITUDE');
    Long=ncread(files(i).name,'LONGITUDE');
    JulDate=ncread(files(i).name,'JULD_LOCATION');
    Depth=[5:5:100,125:25:500,550:50:1000]';
    
    %Remove values with QC flag 4-9
    Temp(Temp_QC>3 & Temp_QC<9)=NaN;
    Sal(Sal_QC>3 & Sal_QC<9)=NaN;        
    Pres(Pres_QC>3 & Pres_QC<9)=NaN;

    %Calculate conserved temperature and absolute salinity using GSW
    %toolbox
    AbsSal=gsw_SA_from_SP(Sal,Pres,Long,Lat);
    ConsTemp=gsw_CT_from_t(AbsSal,Temp,Pres);
    
    %Preallocate arrays
    TempBin=NaN(46,size(Temp,2));
    ConsTempBin=NaN(46,size(ConsTemp,2));
    SalBin=NaN(46,size(Sal,2));
    AbsSalBin=NaN(46,size(AbsSal,2));
    
    %Calculate binned means
    m=size(Lat,1);
    cnt=1;
    for j=1:20
        TempBin(j,:)=mean(Temp(cnt:cnt+4,:),'omitnan');
        ConsTempBin(j,:)=mean(ConsTemp(cnt:cnt+4,:),'omitnan');
        SalBin(j,:)=mean(Sal(cnt:cnt+4,:),'omitnan');
        AbsSalBin(j,:)=mean(AbsSal(cnt:cnt+4,:),'omitnan');
        cnt=cnt+5;
    end
    for j=21:36
        TempBin(j,:)=mean(Temp(cnt:cnt+24,:),'omitnan');
        ConsTempBin(j,:)=mean(ConsTemp(cnt:cnt+24,:),'omitnan');
        SalBin(j,:)=mean(Sal(cnt:cnt+24,:),'omitnan');
        AbsSalBin(j,:)=mean(AbsSal(cnt:cnt+24,:),'omitnan');
        cnt=cnt+25;
    end
    for j=37:46
        TempBin(j,:)=mean(Temp(cnt:cnt+49,:),'omitnan');        
        ConsTempBin(j,:)=mean(ConsTemp(cnt:cnt+49,:),'omitnan');
        SalBin(j,:)=mean(Sal(cnt:cnt+49,:),'omitnan');
        AbsSalBin(j,:)=mean(AbsSal(cnt:cnt+49,:),'omitnan');
        cnt=cnt+50;
    end

    %Create data structure and rearrange data so each row represents one
    %set of measurements (T, S, Dens, etc at a single depth) paired with
    %time, location, TagID, and PTT ID.
    CTD=NaN(46*m,10);
    CTD=array2table(CTD,'VariableNames',{'JulDate','Lat','Long','ID','PTT',...
        'Depth','TempBin','ConsTempBin','SalBin','AbsSalBin'});
    count=1;
    for j=1:m
        CTD.JulDate(count:count+45)=JulDate(j,1);
        CTD.Lat(count:count+45)=Lat(j,1);
        CTD.Long(count:count+45)=Long(j,1);
        CTD.Depth(count:count+45)=Depth;
        CTD.TempBin(count:count+45)=TempBin(:,j)';
        CTD.ConsTempBin(count:count+45)=ConsTempBin(:,j)';
        CTD.SalBin(count:count+45)=SalBin(:,j)';        
        CTD.AbsSalBin(count:count+45)=AbsSalBin(:,j)';
        CTD.ID(count:count+45)=ID;
        CTD.PTT(count:count+45)=ID2;
        count=count+46;
    end
    
    %Concatenate each CTD cast into a single large array.
    CTD_Data_Bin=[CTD_Data_Bin;CTD];  
    clear TempBin Depth Lat Long JulDate CTD cnt SalBin ConsTempBin Pres_QC...
        Temp_QC Sal_QC SFlag TFlag PFlag ID ID2
end

%JulDate in MEOP data is days since 1950-01-01 00:00:00 UTC, add 712237 to
%get MATLAB datenum
CTD_Data_Bin.JulDate=CTD_Data_Bin.JulDate+712237;
CTD_Data_Bin.Year=str2num(datestr(CTD_Data_Bin.JulDate,'yyyy'));
CTD_Data_Bin.Month=str2num(datestr(CTD_Data_Bin.JulDate,'mm'));
CTD_Data_Bin.Long360=CTD_Data_Bin.Long;
CTD_Data_Bin.Long360(CTD_Data_Bin.Long<0)=CTD_Data_Bin.Long(CTD_Data_Bin.Long<0)+360;

%Remove all data from flagged instruments (have bad salinity data)
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==12227 & CTD_Data_Bin.Year==2017)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==12953 & CTD_Data_Bin.Year==2017)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13233)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13235 & CTD_Data_Bin.Year==2016)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13264)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13267 & CTD_Data_Bin.Year==2016)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13267 & CTD_Data_Bin.Year==2017 & CTD_Data_Bin.Month==1)=NaN;
CTD_Data_Bin.SalBin(CTD_Data_Bin.ID==13268 & CTD_Data_Bin.Year==2017)=NaN;

CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==12227 & CTD_Data_Bin.Year==2017)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==12953 & CTD_Data_Bin.Year==2017)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13233)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13235 & CTD_Data_Bin.Year==2016)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13264)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13267 & CTD_Data_Bin.Year==2016)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13267 & CTD_Data_Bin.Year==2017 & CTD_Data_Bin.Month==1)=NaN;
CTD_Data_Bin.AbsSalBin(CTD_Data_Bin.ID==13268 & CTD_Data_Bin.Year==2017)=NaN;
%% Step 2: Load Climatology Data
load('D:\Dropbox\MATLAB\Chapter 1\Data\Climatology\Climatology.mat');

%% Step 3: Match to climatology data
Idx(:,1)=knnsearch(CLon,CTD_Data_Bin.Long);
Idx(:,2)=knnsearch(CLat,CTD_Data_Bin.Lat);
Idx(:,3)=knnsearch(CDepth,CTD_Data_Bin.Depth);
Idx(:,4)=CTD_Data_Bin.Month;

%Preallocate arrays
Temp=NaN(size(CTD_Data_Bin,1),1);
Temp2=NaN(size(CTD_Data_Bin,1),1);
Temp_sd=NaN(size(CTD_Data_Bin,1),1);
Sal=NaN(size(CTD_Data_Bin,1),1);
Sal2=NaN(size(CTD_Data_Bin,1),1);
Sal_sd=NaN(size(CTD_Data_Bin,1),1);

%Pull in Climatology data using indices from knnsearch
for i=1:size(CTD_Data_Bin,1)
    Temp(i)=CTemp(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
    Temp2(i)=CTemp2(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
    Temp_sd(i)=CTsd(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
    Sal(i)=CSal(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
    Sal2(i)=CSal2(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
    Sal_sd(i)=CSsd(Idx(i,1),Idx(i,2),Idx(i,3),Idx(i,4));
end

CTD_Data_Bin.CTemp=Temp;
CTD_Data_Bin.CTemp2=Temp2;
CTD_Data_Bin.CTsd=Temp_sd;
CTD_Data_Bin.CTsd(CTD_Data_Bin.CTsd==0)=NaN;

CTD_Data_Bin.CSal=Sal;
CTD_Data_Bin.CSal2=Sal2;
CTD_Data_Bin.CSsd=Sal_sd;
CTD_Data_Bin.CSsd(CTD_Data_Bin.CSsd==0)=NaN;

%Calculate absolute salinity for climatology
CTD_Data_Bin.CAbsSal=gsw_SA_from_SP(CTD_Data_Bin.CSal,CTD_Data_Bin.Depth,CTD_Data_Bin.Long,CTD_Data_Bin.Lat);
CTD_Data_Bin.CAbsSal2=gsw_SA_from_SP(CTD_Data_Bin.CSal2,CTD_Data_Bin.Depth,CTD_Data_Bin.Long,CTD_Data_Bin.Lat);

%Calculate salinity anomalies
CTD_Data_Bin.AnomS=CTD_Data_Bin.SalBin-CTD_Data_Bin.CSal;
CTD_Data_Bin.AnomS2=CTD_Data_Bin.SalBin-CTD_Data_Bin.CSal2;
CTD_Data_Bin.AnomAbsSal=CTD_Data_Bin.AbsSalBin-CTD_Data_Bin.CAbsSal;
CTD_Data_Bin.AnomAbsSal2=CTD_Data_Bin.AbsSalBin-CTD_Data_Bin.CAbsSal2;

%Calculate conservative temperature for climatology
CTD_Data_Bin.CConsTemp=gsw_CT_from_t(CTD_Data_Bin.CAbsSal,CTD_Data_Bin.CTemp,CTD_Data_Bin.Depth);
CTD_Data_Bin.CConsTemp2=gsw_CT_from_t(CTD_Data_Bin.CAbsSal2,CTD_Data_Bin.CTemp2,CTD_Data_Bin.Depth);

%Calculate temperature anomalies
CTD_Data_Bin.AnomT=CTD_Data_Bin.TempBin-CTD_Data_Bin.CTemp;
CTD_Data_Bin.AnomConsT=CTD_Data_Bin.ConsTempBin-CTD_Data_Bin.CConsTemp;
CTD_Data_Bin.AnomT2=CTD_Data_Bin.ConsTempBin-CTD_Data_Bin.CTemp2;
CTD_Data_Bin.AnomConsT2=CTD_Data_Bin.ConsTempBin-CTD_Data_Bin.CConsTemp2;


%% Step 4: Remove values more or less than 8 standard deviations from the mean 
%from T, S Anomalies (and S T observed values).  Calculate and filter by
%depth layers: Top (<=100m), Mid (101-500m), Bot (501+m)

Numsds=8;

%Calculate sea water density (kg/m^3) using GSW Toolbox
CTD_Data_Bin.Density=gsw_sigma0(CTD_Data_Bin.AbsSalBin,CTD_Data_Bin.ConsTempBin);
CTD_Data_Bin.DeltaD(:)=NaN;
CTD_Data_Bin.DeltaD(2:end)=diff(CTD_Data_Bin.Density);
CTD_Data_Bin.CDens=gsw_sigma0(CTD_Data_Bin.CAbsSal,CTD_Data_Bin.CConsTemp);
CTD_Data_Bin.AnomD=CTD_Data_Bin.Density-CTD_Data_Bin.CDens;

%Separate into depth subsets
Data1=CTD_Data_Bin(CTD_Data_Bin.Depth<=100,:);
Data2=CTD_Data_Bin(CTD_Data_Bin.Depth>100 & CTD_Data_Bin.Depth<=500,:);
Data3=CTD_Data_Bin(CTD_Data_Bin.Depth>500,:);

%Calculate Mean and StDev of T and S for each depth subset.  Multiple StDev
%by Numsds to use as outlier filter
TsdTop=std(Data1.AnomConsT,'omitnan')*Numsds;
SsdTop=std(Data1.AnomS,'omitnan')*Numsds;
TbarTop=mean(Data1.AnomConsT,'omitnan');
SbarTop=mean(Data1.AnomS,'omitnan');

TsdMid=std(Data2.AnomConsT,'omitnan')*Numsds;
SsdMid=std(Data2.AnomS,'omitnan')*Numsds;
TbarMid=mean(Data2.AnomConsT,'omitnan');
SbarMid=mean(Data2.AnomS,'omitnan');

TsdBot=std(Data3.AnomConsT,'omitnan')*Numsds;
SsdBot=std(Data3.AnomS,'omitnan')*Numsds;
TbarBot=mean(Data3.AnomConsT,'omitnan');
SbarBot=mean(Data3.AnomS,'omitnan');

%Remove measurements more than +/- Numsds from means
Data1(Data1.AnomConsT>(TbarTop+TsdTop),:)=[];
Data1(Data1.AnomConsT<(TbarTop-TsdTop),:)=[];
Data2(Data2.AnomConsT>(TbarMid+TsdMid),:)=[];
Data2(Data2.AnomConsT<(TbarMid-TsdMid),:)=[];
Data3(Data3.AnomConsT>(TbarBot+TsdBot),:)=[];
Data3(Data3.AnomConsT<(TbarBot-TsdBot),:)=[];

Data1(Data1.AnomS>(SbarTop+SsdTop),:)=[];
Data1(Data1.AnomS<(SbarTop-SsdTop),:)=[];
Data2(Data2.AnomS>(SbarMid+SsdMid),:)=[];
Data2(Data2.AnomS<(SbarMid-SsdMid),:)=[];
Data3(Data3.AnomS>(SbarBot+SsdBot),:)=[];
Data3(Data3.AnomS<(SbarBot-SsdBot),:)=[];

%Recombine filtered data
CTD_Data_Bin2=[Data1;Data2;Data3];

%post-filter histograms  
fig=figure(1);
fig.Color = [1 1 1];
t=tiledlayout(2,2);
t.TileSpacing='compact';
nexttile
ax=gca;
histogram(CTD_Data_Bin.AnomT);
ax.YLim=[0 50];
ax.Title.String='Temperature Anomalies';
nexttile
ax=gca;
histogram(CTD_Data_Bin.AnomS);
ax.YLim=[0 50];
ax.Title.String='Salinity Anomalies';
nexttile
ax=gca;
histogram(CTD_Data_Bin.AnomD);
ax.YLim=[0 50];
ax.Title.String='Density Anomalies';

clear Data1 Data2 Data3

%% Step 5: Calculate sigmaTa for each location

SigmaTagT=0.02; %from Roquet et al. 2011
CTD_Data_Bin.SigmaAnomT = sqrt(SigmaTagT^2 + CTD_Data_Bin.CTsd.^2);

SigmaTagS=0.03; %from Roquet et al. 2011
CTD_Data_Bin.SigmaAnomS = sqrt(SigmaTagS^2 + CTD_Data_Bin.CSsd.^2);

CTD_Data_Bin.AnomTSigmas=CTD_Data_Bin.AnomT./CTD_Data_Bin.SigmaAnomT;
CTD_Data_Bin.AnomConsTSigmas=CTD_Data_Bin.AnomConsT./CTD_Data_Bin.SigmaAnomT;
CTD_Data_Bin.AnomSSigmas=CTD_Data_Bin.AnomS./CTD_Data_Bin.SigmaAnomS;

CTD_Data_Bin.AnomT2Sigmas=CTD_Data_Bin.AnomT2./CTD_Data_Bin.SigmaAnomT;
CTD_Data_Bin.AnomConsT2Sigmas=CTD_Data_Bin.AnomConsT2./CTD_Data_Bin.SigmaAnomT;
CTD_Data_Bin.AnomS2Sigmas=CTD_Data_Bin.AnomS2./CTD_Data_Bin.SigmaAnomS;

%% Step 6: Calculate spiciness using GSW toolbox

CTD_Data_Bin.Spice=gsw_spiciness0(CTD_Data_Bin.AbsSalBin,CTD_Data_Bin.ConsTempBin);
CTD_Data_Bin.CSpice=gsw_spiciness0(CTD_Data_Bin.CAbsSal,CTD_Data_Bin.CConsTemp);
CTD_Data_Bin.SpiceAnom=CTD_Data_Bin.Spice-CTD_Data_Bin.CSpice;

%% Step 7: Save data as .mat file for analysis
save('All_CTD_Data_Bin_V2.mat','-v7.3','CTD_Data_Bin');
clear Idx Temp count files Temp_sd Sal Sal_sd
