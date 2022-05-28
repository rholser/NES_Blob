%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 25-May-2022

%%%%%Required scripts, functions, and toolboxes:
    %GSW toolbox (TEOS 10)
    %WOA18 climatology netCDF files
    %All_CTD_Data_Bin_V2.mat (from Import_CalcAnomalyBinned_MEOP_V2.m)
%%%%%%%

%Finds the latitude and longitude of nearest climatology data that matches
%the T S signature in seal data for each grid cell of a given density layer
%by bootstrapping a cost function (see Edwards et al. 2015).  
%Saves the results as .mat file for easy recall.  Bootstrapping step can
%take 12+ hours depending on number of iterations and computing power.

%Creates plots of:
    %Spice anomaly contours with arrows
    %Depth anomaly contours
    %Climatology values
    
%Plots generated are not automatically saved - must be exported manually.

%% Pre-define density surface and lat/lon range to analyze
clear

dens=27;
int=0.075; %interval around target density to include
LatMin=20;
LatMax=63;
LonMin=165;
LonMax=245;
LatInt=5;
LonInt=5;

%Define Grid for seal data
LatMin2=28.5;
LatMax2=62.5;
LonMin2=178.5;
LonMax2=242.5;

%% Step 1: Prepare Climatoloty Data

%Load Annual WOA18 Temp and Sal on 1x1 grid
files=dir('*.nc');

CLat=double(ncread(files(1).name,'lat'));
CLon=double(ncread(files(1).name,'lon'));
CDepth=ncread(files(1).name,'depth');
CSal(:,:,:)=ncread(files(1).name,'s_an');
CSal2(:,:,:)=ncread(files(1).name,'s_mn');
CSsd(:,:,:)=ncread(files(1).name,'s_sd');
CTemp(:,:,:)=ncread(files(3).name,'t_an');
CTemp2(:,:,:)=ncread(files(3).name,'t_mn');
CTsd(:,:,:)=ncread(files(3).name,'t_sd');
CDepthInt=(0:5:5500)';

%Create 0-360 longitude structure
CLon2=CLon;
CLon2(CLon2<0)=CLon2(CLon2<0)+360;
%Longitude index for climatology data
LonInd(:,1)=find(CLon2>=LonMin & CLon2<=LonMax);
LonInd(:,2)=CLon2(CLon2>=LonMin & CLon2<=LonMax); 
LonInd=sortrows(LonInd,2);
%Latitude index for climatology data
LatInd(:,1)=find(CLat>=LatMin & CLat<=LatMax);
LatInd(:,2)=CLat(CLat>=LatMin & CLat<=LatMax);
%Lat/lon index for seal data
LatInd2=LatInd(LatInd(:,2)>=LatMin2 & LatInd(:,2)<=LatMax2,:);
LonInd2=LonInd(LonInd(:,2)>=LonMin2 & LonInd(:,2)<=LonMax2,:);

CSalInt=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
CTempInt=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
CAbsSal=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
CConsTemp=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));

CSsdInt=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
CTsdInt=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
CDensInt=NaN(size(LonInd,1),size(LatInd,1),size(CDepthInt,1));
DTemp=NaN(size(LonInd,1),size(LatInd,1));
DSal=NaN(size(LonInd,1),size(LatInd,1));
DDepth=NaN(size(LonInd,1),size(LatInd,1));
DDens=NaN(size(LonInd,1),size(LatInd,1));
for i=1:size(LonInd,1)
    for j=1:size(LatInd,1)
        CSalInt(i,j,:)=interp1(CDepth,squeeze(CSal(LonInd(i,1),LatInd(j,1),:)),CDepthInt);
        CTempInt(i,j,:)=interp1(CDepth,squeeze(CTemp(LonInd(i,1),LatInd(j,1),:)),CDepthInt);
        CSsdInt(i,j,:)=interp1(CDepth,squeeze(CSsd(LonInd(i,1),LatInd(j,1),:)),CDepthInt);
        CTsdInt(i,j,:)=interp1(CDepth,squeeze(CTsd(LonInd(i,1),LatInd(j,1),:)),CDepthInt);

        CConsTemp(i,j,:)=gsw_CT_from_t(squeeze(CSalInt(i,j,:)),squeeze(CTempInt(i,j,:)),CDepthInt(:,1));
        CAbsSal(i,j,:)=gsw_SA_from_SP(squeeze(CSalInt(i,j,:)),CDepthInt(:,1),LonInd(j,2),LatInd(j,2));
        CDensInt(i,j,:)=gsw_sigma0(squeeze(CAbsSal(i,j,:)),squeeze(CConsTemp(i,j,:)));

        [~,ind]=min(abs(CDensInt(i,j,:)-dens));
        %DTemp(i,j)=CTempInt(i,j,ind);
        %DSal(i,j)=CSalInt(i,j,ind);
        DTemp(i,j)=CConsTemp(i,j,ind);
        DSal(i,j)=CAbsSal(i,j,ind);
        DDepth(i,j)=CDepthInt(ind,1);
        DDens(i,j)=CDensInt(i,j,ind);
        DTsd(i,j)=CTsdInt(i,j,ind);
        DSsd(i,j)=CSsdInt(i,j,ind);
    end
end

%Reshape data and narrow to desired lat/lon range
DTemp=DTemp';
DDens=DDens';
DDepth=DDepth';
DSal=DSal';
DTsd=DTsd';
DSsd=DSsd';

clear CLon CLon2 CLat CSal CSal2 CTemp CTemp2 CDepth CTsd CSsd CDens CDens2 files i j
%% Step 2: Prepare Eseal Data

%Load data and subset to density
load('All_CTD_Data_Bin_V2.mat');
Data=CTD_Data_Bin(CTD_Data_Bin.Density>=(dens-int) & CTD_Data_Bin.Density<=(dens+int),:);

%Subset by years
Data14=Data(Data.Year==2014,:);
Data15=Data(Data.Year==2015,:);
Data16=Data(Data.Year==2016,:);
Data17=Data(Data.Year==2017,:);

%Calculate mean values of T and S for each year in each grid cell
TGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
TGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
TGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
TGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));
SGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
SGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
SGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
SGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));
DGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
DGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
DGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
DGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));

TaGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
TaGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
TaGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
TaGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));
SaGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
SaGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
SaGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
SaGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));
SpGrid_14=NaN(size(LatInd2,1),size(LonInd2,1));
SpGrid_15=NaN(size(LatInd2,1),size(LonInd2,1));
SpGrid_16=NaN(size(LatInd2,1),size(LonInd2,1));
SpGrid_17=NaN(size(LatInd2,1),size(LonInd2,1));

for i=1:size(LatInd2,1)-1
    temp14=Data14(Data14.Lat>=LatInd2(i,2) & Data14.Lat<LatInd2(i+1,2),:);
    temp15=Data15(Data15.Lat>=LatInd2(i,2) & Data15.Lat<LatInd2(i+1,2),:);
    temp16=Data16(Data16.Lat>=LatInd2(i,2) & Data16.Lat<LatInd2(i+1,2),:);
    temp17=Data17(Data17.Lat>=LatInd2(i,2) & Data17.Lat<LatInd2(i+1,2),:);
    for j=1:size(LonInd2,1)-1
        temp14b=temp14(temp14.Long360>=LonInd2(j,2) & temp14.Long360<LonInd2(j+1,2),:);
        temp15b=temp15(temp15.Long360>=LonInd2(j,2) & temp15.Long360<LonInd2(j+1,2),:);
        temp16b=temp16(temp16.Long360>=LonInd2(j,2) & temp16.Long360<LonInd2(j+1,2),:);
        temp17b=temp17(temp17.Long360>=LonInd2(j,2) & temp17.Long360<LonInd2(j+1,2),:);
        
        TGrid_14(i,j)=mean(temp14b.ConsTemp,'omitnan');
        SGrid_14(i,j)=mean(temp14b.AbsSal,'omitnan');
        TGrid_15(i,j)=mean(temp15b.ConsTemp,'omitnan');
        SGrid_15(i,j)=mean(temp15b.AbsSal,'omitnan');
        TGrid_16(i,j)=mean(temp16b.ConsTemp,'omitnan');
        SGrid_16(i,j)=mean(temp16b.AbsSal,'omitnan');
        TGrid_17(i,j)=mean(temp17b.ConsTemp,'omitnan');
        SGrid_17(i,j)=mean(temp17b.AbsSal,'omitnan');
        DGrid_14(i,j)=mean(temp14b.Depth,'omitnan');
        DGrid_15(i,j)=mean(temp15b.Depth,'omitnan');
        DGrid_16(i,j)=mean(temp16b.Depth,'omitnan');
        DGrid_17(i,j)=mean(temp17b.Depth,'omitnan');
        
        SpGrid_14(i,j)=mean(temp14b.SpiceAnom,'omitnan');
        SpGrid_15(i,j)=mean(temp15b.SpiceAnom,'omitnan');
        SpGrid_16(i,j)=mean(temp16b.SpiceAnom,'omitnan');
        SpGrid_17(i,j)=mean(temp17b.SpiceAnom,'omitnan');
        
        TaGrid_14(i,j)=mean(temp14b.AnomConsT,'omitnan');
        SaGrid_14(i,j)=mean(temp14b.AnomAbsSal,'omitnan');
        TaGrid_15(i,j)=mean(temp15b.AnomConsT,'omitnan');
        SaGrid_15(i,j)=mean(temp15b.AnomAbsSal,'omitnan');
        TaGrid_16(i,j)=mean(temp16b.AnomConsT,'omitnan');
        SaGrid_16(i,j)=mean(temp16b.AnomAbsSal,'omitnan');
        TaGrid_17(i,j)=mean(temp17b.AnomConsT,'omitnan');
        SaGrid_17(i,j)=mean(temp17b.AnomAbsSal,'omitnan');
        clear temp14b temp15b temp16b temp17b
    end
    clear temp14 temp15 temp16 temp17
end


SpGrid_Diff_15=SpGrid_15-SpGrid_14;
SpGrid_Diff_16=SpGrid_16-SpGrid_15;
SpGrid_Diff_17=SpGrid_17-SpGrid_16;

%Calculate T sigma and S sigma based on presence of seal data (used 2014 as
%broadest coverage)
DSsd=DSsd(LatInd(:,2)>=LatMin2 & LatInd(:,2)<=LatMax2,LonInd(:,2)>=LonMin2 & LonInd(:,2)<=LonMax2);
DTsd=DTsd(LatInd(:,2)>=LatMin2 & LatInd(:,2)<=LatMax2,LonInd(:,2)>=LonMin2 & LonInd(:,2)<=LonMax2);

for i=1:size(LatInd2,1)
    for j=1:size(LonInd2,1)
        if isnan(TGrid_14(i,j))
            DTsd(i,j)=NaN;
            DSsd(i,j)=NaN;
        end
    end
end

DSsd_Mean=mean(DSsd,'all','omitnan');
DTsd_Mean=mean(DTsd,'all','omitnan');
DSsd_sd=std(DSsd,0,'all','omitnan');
DTsd_sd=std(DTsd,0,'all','omitnan');

clear Data Data14 Data15 Data16 Data17

%% Step 3: Calculate Depth Anomaly
%rearrange data for continuous 360 longitudes
DDepth2=DDepth(9:end,14:78);

DaGrid_14=DGrid_14-DDepth2;
DaGrid_15=DGrid_15-DDepth2;
DaGrid_16=DGrid_16-DDepth2;
DaGrid_17=DGrid_17-DDepth2;

%% Step 4: Calculate Closest Climatology Cell Using Cost Function

%Create arrays to store final lat-long coordinates (lat x long x 2 array).  
%Layer 1 latitude, layer 2 longitude of closest matching point on
%climatology grid.

N=1000; %Number of bootstrap iterations

Dist_14=NaN(size(LatInd2,1),size(LonInd2,1),3,N);
Dist_15=NaN(size(LatInd2,1),size(LonInd2,1),3,N);
Dist_16=NaN(size(LatInd2,1),size(LonInd2,1),3,N);
Dist_17=NaN(size(LatInd2,1),size(LonInd2,1),3,N);

s=rng('default');

%Generate N random, normally distributed perterbations around mean of 0
%using the mean stadnard deviations calculated above
Trand=DTsd_Mean*randn(N,1);
Srand=DSsd_Mean*randn(N,1);
    
%perturb both T and S in each iteration and calculate closest match within 
%10x10 degree box for each.
tic
for n=1:N
    n
    for i=1:size(LatInd2,1)-1
        for j=1:size(LonInd2,1)-1
            temp14=NaN(size(DTemp,1),size(DTemp,2));
            temp15=NaN(size(DTemp,1),size(DTemp,2));
            temp16=NaN(size(DTemp,1),size(DTemp,2));
            temp17=NaN(size(DTemp,1),size(DTemp,2));
            
            T14=TGrid_14(i,j)+Trand(n,1);
            S14=SGrid_14(i,j)+Srand(n,1);
            T15=TGrid_15(i,j)+Trand(n,1);
            S15=SGrid_15(i,j)+Srand(n,1);
            T16=TGrid_16(i,j)+Trand(n,1);
            S16=SGrid_16(i,j)+Srand(n,1);
            T17=TGrid_17(i,j)+Trand(n,1);
            S17=SGrid_17(i,j)+Srand(n,1);
            
            for k=1:size(LatInd,1)-1
                for l=1:size(LonInd,1)-1
                    %Calculate Temp and Sal distance between seal values and
                    %climatology, normalized by climatology mn field sigma.
                    %Preturb T and S using Trand and Srand.
                    temp14(k,l)=((T14-DTemp(k,l))/DTsd_Mean)^2 + ...
                        ((S14-DSal(k,l))/DSsd_Mean)^2;               
                    temp15(k,l)=((T15-DTemp(k,l))/DTsd_Mean)^2 + ...
                        ((S15-DSal(k,l))/DSsd_Mean)^2;              
                    temp16(k,l)=((T16-DTemp(k,l))/DTsd_Mean)^2 + ...
                        ((S16-DSal(k,l))/DSsd_Mean)^2;
                    temp17(k,l)=((T17-DTemp(k,l))/DTsd_Mean)^2 + ...
                        ((S17-DSal(k,l))/DSsd_Mean)^2;
                end
            end

            iLatMin=find(LatInd(:,2)==(LatInd2(i,2)-LatInt));
                if isempty(iLatMin)
                    iLatMin=1;
                end
            iLatMax=find(LatInd(:,2)==(LatInd2(i,2)+LatInt));
                if isempty(iLatMax)
                    iLatMax=size(LatInd,1);
                end
            iLonMin=find(LonInd(:,2)==(LonInd2(j,2)-LonInt));
                 if isempty(iLonMin)
                    iLonMin=1;
                end
            iLonMax=find(LonInd(:,2)==(LonInd2(j,2)+LonInt));
                if isempty(iLonMax)
                    iLonMax=size(LonInd,1);
                end

            MinVal_14=min(temp14(iLatMin:iLatMax,iLonMin:iLonMax),[],'all');
            if ~isnan(MinVal_14)
                [row14,col14]=find(temp14(iLatMin:iLatMax,iLonMin:iLonMax)==MinVal_14);
                Dist_14(i,j,1,n)=LatInd(row14+iLatMin-1,2);
                Dist_14(i,j,2,n)=LonInd(col14+iLonMin-1,2);
                Dist_14(i,j,3,n)=MinVal_14;
            end
            MinVal_15=min(temp15(iLatMin:iLatMax,iLonMin:iLonMax),[],'all');
            if ~isnan(MinVal_15)
                [row15,col15]=find(temp15(iLatMin:iLatMax,iLonMin:iLonMax)==MinVal_15);
                Dist_15(i,j,1,n)=LatInd(row15+iLatMin-1,2);
                Dist_15(i,j,2,n)=LonInd(col15+iLonMin-1,2);
                Dist_15(i,j,3,n)=MinVal_15;
            end

            MinVal_16=min(temp16(iLatMin:iLatMax,iLonMin:iLonMax),[],'all');
            if ~isnan(MinVal_16)
                [row16,col16]=find(temp16(iLatMin:iLatMax,iLonMin:iLonMax)==MinVal_16);
                Dist_16(i,j,1,n)=LatInd(row16+iLatMin-1,2);
                Dist_16(i,j,2,n)=LonInd(col16+iLonMin-1,2);
                Dist_16(i,j,3,n)=MinVal_16;
            end

            MinVal_17=min(temp17(iLatMin:iLatMax,iLonMin:iLonMax),[],'all');
            if ~isnan(MinVal_17)
                [row17,col17]=find(temp17(iLatMin:iLatMax,iLonMin:iLonMax)==MinVal_17);
                Dist_17(i,j,1,n)=LatInd(row17+iLatMin-1,2);
                Dist_17(i,j,2,n)=LonInd(col17+iLonMin-1,2);
                Dist_17(i,j,3,n)=MinVal_17;
            end
            clear row14 col14 MinVal_14 row15 col15 MinVal_15 row16 col16...
                MinVal_16 row17 col17 MinVal_17 temp14 temp15 temp16 temp17...
                T14 T15 T16 T17 S14 S15 S16 S17
        end
    end
end
toc

Dist_14_Mean=mean(Dist_14,4,'omitnan');
Dist_15_Mean=mean(Dist_15,4,'omitnan');
Dist_16_Mean=mean(Dist_16,4,'omitnan');
Dist_17_Mean=mean(Dist_17,4,'omitnan');

Dist_14_Median=median(Dist_14,4,'omitnan');
Dist_15_Median=median(Dist_15,4,'omitnan');
Dist_16_Median=median(Dist_16,4,'omitnan');
Dist_17_Median=median(Dist_17,4,'omitnan');

%save resulting data arrays - file naming is not automated but should reflect 
%isopycnal and number of iterations used
save('DensitySurfaceDistance_27_N1000_V3.mat','-v7.3',...
    'Dist_14','Dist_15','Dist_16','Dist_17',...
    'Dist_14_Mean','Dist_15_Mean','Dist_16_Mean','Dist_17_Mean');

%% Step 5: Plots

%Spice figures with Arrows
R = georasterref('RasterSize', size(SpGrid_14(:,:)),'Latlim', [min(LatInd2(:,2)) max(LatInd2(:,2))],...
    'Lonlim', [min(LonInd2(:,2)) max(LonInd2(:,2))]);

%define spice scale based on density being displayed
if dens==27
    spmax=0.3;
    spmin=-0.3;
    spint=0.05;
else
    spmax=0.8;
    spmin=-0.8;
    spint=0.1;
end

figure; %2014
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',40);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(SpGrid_14,R,[-1:spint:1],'LineStyle','none');
ax.Colormap=([customcolormap_preset('pasteljet')]);
caxis([spmin,spmax]);
c=contourcbar;
c.Location='southoutside';
c.FontSize=40;
c.Label.String = 'Spice Anomaly';
c.Color='k';
c.LineWidth=1;
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
for i=1:size(LatInd2,1)
    for j=1:2:size(LonInd2,1)
        [XX,YY] = arrow([Dist_14_Mean(i,j,1);LatInd2(i,2)],[Dist_14_Mean(i,j,2);LonInd2(j,2)]);
        plotm([[LatInd2(i,2);Dist_14_Mean(i,j,1)] ; NaN*ones(1,size([LatInd2(i,2);Dist_14_Mean(i,j,1)],2)) ; XX], ...
          [[LonInd2(j,2); Dist_14_Mean(i,j,2)] ; NaN*ones(1,size([LonInd2(j,2); Dist_14_Mean(i,j,2)],2)) ; YY],'k')
    end
end

figure; %2015
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',40);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(SpGrid_Diff_15,R,[-1:spint:1],'LineStyle','none');
ax.Colormap=([customcolormap_preset('pasteljet')]);
caxis([spmin,spmax]);
c=contourcbar;
c.Location='southoutside';
c.FontSize=40;
c.Label.String = 'Spice Anomaly';
c.Color='k';
c.LineWidth=1;
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
for i=1:size(LatInd2,1)
    for j=1:2:size(LonInd2,1)
        [XX,YY] = arrow([Dist_15_Mean(i,j,1);LatInd2(i,2)],[Dist_15_Mean(i,j,2);LonInd2(j,2)]);
        plotm([[LatInd2(i,2);Dist_15_Mean(i,j,1)] ; NaN*ones(1,size([LatInd2(i,2);Dist_15_Mean(i,j,1)],2)) ; XX], ...
          [[LonInd2(j,2); Dist_15_Mean(i,j,2)] ; NaN*ones(1,size([LonInd2(j,2); Dist_15_Mean(i,j,2)],2)) ; YY],'k')
    end
end

figure; %2016
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',40);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(SpGrid_16,R,[-1:spint:1],'LineStyle','none');
ax.Colormap=([customcolormap_preset('pasteljet')]);
caxis([spmin,spmax]);
c=contourcbar;
c.Location='southoutside';
c.FontSize=40;
c.Label.String = 'Spice Anomaly';
c.Color='k';
c.LineWidth=1;
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
for i=1:size(LatInd2,1)
    for j=1:2:size(LonInd2,1)
        [XX,YY] = arrow([Dist_16_Mean(i,j,1);LatInd2(i,2)],[Dist_16_Mean(i,j,2);LonInd2(j,2)]);
        plotm([[LatInd2(i,2);Dist_16_Mean(i,j,1)] ; NaN*ones(1,size([LatInd2(i,2);Dist_16_Mean(i,j,1)],2)) ; XX], ...
          [[LonInd2(j,2); Dist_16_Mean(i,j,2)] ; NaN*ones(1,size([LonInd2(j,2); Dist_16_Mean(i,j,2)],2)) ; YY],'k')
    end
end

figure; %2017
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',40);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(SpGrid_17,R,[-1:spint:1],'LineStyle','none');
ax.Colormap=([customcolormap_preset('pasteljet')]);
caxis([spmin,spmax]);
c=contourcbar;
c.Location='southoutside';
c.FontSize=40;
c.Label.String = 'Spice Anomaly';
c.Color='k';
c.LineWidth=1;
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
for i=1:size(LatInd2,1)
    for j=1:2:size(LonInd2,1)
        [XX,YY] = arrow([Dist_17_Mean(i,j,1);LatInd2(i,2)],[Dist_17_Mean(i,j,2);LonInd2(j,2)]);
        plotm([[LatInd2(i,2);Dist_17_Mean(i,j,1)] ; NaN*ones(1,size([LatInd2(i,2);Dist_17_Mean(i,j,1)],2)) ; XX], ...
          [[LonInd2(j,2); Dist_17_Mean(i,j,2)] ; NaN*ones(1,size([LonInd2(j,2); Dist_17_Mean(i,j,2)],2)) ; YY],'k')
    end
end

%% Depth Anomaly Plots
R = georasterref('RasterSize', size(TGrid_14(:,:)),'Latlim', [min(LatInd2(:,2)) max(LatInd2(:,2))],...
    'Lonlim', [min(LonInd2(:,2)) max(LonInd2(:,2))]);

figure; %2014
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DaGrid_14,R,[-100:10:100],'LineStyle','none');
caxis([-100,100]);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
ax.Colormap=flipud(customcolormap_preset('red-yellow-blue'));
c=contourcbar;
c.Location='southoutside';
c.FontSize=36;
c.Label.String = 'Depth Anomaly (m)';
c.Color='k';
c.LineWidth=1;

figure; %2015
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.Colormap=flipud(parula);
ax.YColor='w';
ax.XColor='w';
contourfm(DaGrid_15,R,[-100:10:100],'LineStyle','none');
caxis([-100,100]);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
ax.Colormap=flipud(customcolormap_preset('red-yellow-blue'));
c=contourcbar;
c.Location='southoutside';
c.FontSize=36;
c.Label.String = 'Depth Anomaly (m)';
c.Color='k';
c.LineWidth=1;

figure; %2016
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DaGrid_16,R,[-100:10:100],'LineStyle','none');
caxis([-100,100]);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
ax.Colormap=flipud(customcolormap_preset('red-yellow-blue'));
c=contourcbar;
c.Location='southoutside';
c.FontSize=36;
c.Label.String = 'Depth Anomaly (m)';
c.Color='k';
c.LineWidth=1;

figure; %2017
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[178 242],'maplatlimit',[28,61],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DaGrid_17,R,[-100:10:100],'LineStyle','none');
caxis([-100,100]);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
ax.Colormap=flipud(customcolormap_preset('red-yellow-blue'));
c=contourcbar;
c.Location='southoutside';
c.FontSize=36;
c.Label.String = 'Depth Anomaly (m)';
c.Color='k';
c.LineWidth=1;

%% Climatology Plots

%Plot to check data properly subset
R = georasterref('RasterSize', size(DSal),'Latlim', [min(LatInd(:,2)) max(LatInd(:,2))],...
    'Lonlim', [min(LonInd(:,2)) max(LonInd(:,2))]);

figure;
fig.Color = [1 1 1];
axesm('robinson','maplonlimit',[min(LonInd(:,2)) max(LonInd(:,2))],'maplatlimit',[min(LatInd(:,2)) max(LatInd(:,2))],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
geoshow(DTsd,R,'DisplayType','texturemap')%use 'surface' so that NaNs can be white
%contourfm(DTempsd,R,[0:2:18],'LineStyle','none');
ax.Colormap=parula;
caxis([0,2.5]);
c=contourcbar;
c.Location='southoutside';
c.FontSize=24;
c.Label.String = 'Temperature \circ C';
c.Color='k';
c.LineWidth=1;
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);

%Climatology Temperature
figure;
axesm('robinson','maplonlimit',[LonMin LonMax],'maplatlimit',[LatMin,LatMax],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DTemp,R,[0:0.5:18],'LineStyle','none');
ax.Colormap=turbo;
c=contourcbar;
caxis([1,11]);
c.Location='southoutside';
c.FontSize=24;
c.Label.String = 'Temperature \circ C';
c.Color='k';
c.LineWidth=1;
[C,h]=contourm(LatInd(:,2),LonInd(:,2),DSal,...
    'k','LevelStep',0.4,'LineWidth',1.5,'LineStyle','--');
t=clabelm(C,h);
[C,h]=contourm(LatInd(:,2),LonInd(:,2),DSal,...
    'k','LevelStep',0.2,'LineWidth',1.5,'LineStyle',':');
[C,h]=contourm(LatInd(:,2),LonInd(:,2),DSal,...
    'k','LevelStep',1,'LineWidth',1.5);
t=clabelm(C,h);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);

%Climatology Salinity
figure;
axesm('robinson','maplonlimit',[LonMin LonMax],'maplatlimit',[LatMin,LatMax],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DSal,R,[32.4:0.2:34.8],'LineStyle','none');
ax.Colormap=parula;
c=contourcbar;
caxis([32.4,34.8]);
c.Location='southoutside';
c.Ticks=[32.4:0.2:34.8]
c.FontSize=24;
c.Label.String = 'Salinity';
c.Color='k';
c.LineWidth=1;
contourm(LatInd(:,2),LonInd(:,2),DTemp,...
    'k','LevelStep',1,'LineWidth',1.5,'LineStyle',':');
[C,h]=contourm(LatInd(:,2),LonInd(:,2),DTemp,...
    'k','LevelStep',2,'LineWidth',1.5);
t=clabelm(C,h);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);

%Climatology Depth
figure;
axesm('robinson','maplonlimit',[LonMin LonMax],'maplatlimit',[LatMin,LatMax],'frame','on',...
    'Grid','on','MeridianLabel','on','ParallelLabel','on','Fontsize',36);
ax=gca;
ax.Box='off';
ax.YColor='w';
ax.XColor='w';
contourfm(DDepth,R,[50:25:600],'LineStyle','none');
ax.Colormap=flipud(parula);
c=contourcbar;
caxis([50,600]);
c.Location='southoutside';
c.FontSize=24;
c.Label.String = 'Depth (m)';
c.Color='k';
c.LineWidth=1;
contourm(LatInd(:,2),LonInd(:,2),DTemp,...
    'k','LevelStep',1,'LineWidth',1.5,'LineStyle',':');
[C,h]=contourm(LatInd(:,2),LonInd(:,2),DTemp,...
    'k','LevelStep',2,'LineWidth',1.5);
t=clabelm(C,h);
geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
