%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 28-May-2022

%%%%% Required scripts, functions, and toolboxes:
    %GSW toolbox (TEOS 10)
    %Climatology.mat file (from Import_Climatology_Data.m)
    %WOA18 annual climatology netCDF files
    %All_CTD_Data_Bin_V2.mat (from Import_CalcAnomalyBinned_MEOP_V2.m)
%%%%%%%

%Script to calculate and plot seasonal cycle in spice anomaly.

%% Step 1: Import monthly climatology data, calculate TCons, AbsSal, and Spice. 

%Set Latitude-longitude bounds for 40-50Nx130-150W
%Use 5m depth data interpolation

Depths=(5:5:1000)';
LatMax=50;
LatMin=40;
LonMin=-150; 
LonMax=-130;
Months=1:48;

load('Climatology.mat');

CAbsSal=NaN(360,180,57,12);
CAbsSal2=NaN(360,180,57,12);
CConsTemp=NaN(360,180,57,12);
CConsTemp2=NaN(360,180,57,12);
CSpice=NaN(360,180,57,12);
CDensity=NaN(360,180,57,12);

for i=1:12
    %Calculate absolute salinity for climatology
    for j=1:size(CSal,1)
        CAbsSal(j,:,:,i)=gsw_SA_from_SP(squeeze(CSal(j,:,:,i)),CDepth,CLon(j),CLat);
        CAbsSal2(j,:,:,i)=gsw_SA_from_SP(squeeze(CSal2(j,:,:,i)),CDepth,CLon(j),CLat);

        %Calculate conservative temperature for climatology
        CConsTemp(j,:,:,i)=gsw_CT_from_t(squeeze(CAbsSal(j,:,:,i)),squeeze(CTemp(j,:,:,i)),CDepth);
        CConsTemp2(j,:,:,i)=gsw_CT_from_t(squeeze(CAbsSal2(j,:,:,i)),squeeze(CTemp2(j,:,:,i)),CDepth);
        
        %Calculate density from absolute salinity and conservative temp
        CDensity=gsw_sigma0(CAbsSal,CConsTemp);

        %Calculate spiciness
        CSpice(j,:,:,i)=gsw_spiciness0(squeeze(CAbsSal(j,:,:,i)),squeeze(CConsTemp(j,:,:,i)));
    end
end

%% Step 2: Import annual climatology data, calculate TCons, AbsSal, and Spice. 

files=dir('*.nc');

CAnLat=double(ncread(files(1).name,'lat'));
CAnLon=double(ncread(files(1).name,'lon'));
CAnDepth=ncread(files(1).name,'depth');
CAnSal(:,:,:)=ncread(files(1).name,'s_an');
CAnSal2(:,:,:)=ncread(files(1).name,'s_mn');
CAnSsd(:,:,:)=ncread(files(1).name,'s_sd');
CAnTemp(:,:,:)=ncread(files(3).name,'t_an');
CAnTemp2(:,:,:)=ncread(files(3).name,'t_mn');
CAnTsd(:,:,:)=ncread(files(3).name,'t_sd');

CAnAbsSal=NaN(360,180,102);
CAnAbsSal2=NaN(360,180,102);
CAnConsTemp=NaN(360,180,102);
CAnConsTemp2=NaN(360,180,102);
CAnSpice=NaN(360,180,102);

for j=1:size(CSal,1)
    CAnAbsSal(j,:,:)=gsw_SA_from_SP(squeeze(CAnSal(j,:,:)),CAnDepth,CAnLon(j),CAnLat);
    CAnAbsSal2(j,:,:)=gsw_SA_from_SP(squeeze(CAnSal2(j,:,:)),CAnDepth,CAnLon(j),CAnLat);

    %Calculate conservative temperature for climatology
    CAnConsTemp(j,:,:)=gsw_CT_from_t(squeeze(CAnAbsSal(j,:,:)),squeeze(CAnTemp(j,:,:)),CAnDepth);
    CAnConsTemp2(j,:,:)=gsw_CT_from_t(squeeze(CAnAbsSal2(j,:,:)),squeeze(CAnTemp2(j,:,:)),CAnDepth);

    %Calculate spiciness
    CAnSpice(j,:,:)=gsw_spiciness0(squeeze(CAnAbsSal(j,:,:)),squeeze(CAnConsTemp(j,:,:)));
end

%% Step 3: subtract annual spice from all monthly spice values and calculate
% mean monthly spice within box, by depth

CSpiceAnom=NaN(360,180,57,12);
CSpAnomBox=NaN(57,12);
CDensBox=NaN(57,12);

for i=1:12
   CSpiceAnom(:,:,:,i)=CSpice(:,:,:,i)-CAnSpice(:,:,1:57); 
   for j=1:57
        CDensBox(j,i)=mean(CDensity(find(ceil(CLon)==LonMin):find(ceil(CLon)==LonMax),...
            find(ceil(CLat)==LatMin):find(ceil(CLat)==LatMax),j,i),'all','omitnan');
        CSpAnomBox(j,i)=mean(CSpiceAnom(find(ceil(CLon)==LonMin):find(ceil(CLon)==LonMax),...
            find(ceil(CLat)==LatMin):find(ceil(CLat)==LatMax),j,i),'all','omitnan');
   end
   
end
CSpAnomBox(:,13:24)=CSpAnomBox(:,1:12);
CSpAnomBox(:,25:48)=CSpAnomBox(:,1:24);
CDensBox(:,13:24)=CDensBox(:,1:12);
CDensBox(:,25:48)=CDensBox(:,1:24);

%% Step 4: Interpolate to 5m depth
CSpAnomInt=NaN(200,48);
CDensInt=NaN(200,48);

for i=1:48
    CSpAnomInt(:,i)=interp1(CDepth,CSpAnomBox(:,i),Depths);
    CDensInt(:,i)=interp1(CDepth,CDensBox(:,i),Depths);
end

%% Step 5: Create plot of annual cycle x 4 to match observed data

fig=figure(1);
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Months,Depths,CSpAnomInt,[-0.6 -.4 -0.2 -0.1 -0.05 0 0.05 0.1 0.2 0.4 0.6]);
hold on
ax.CLim=[-0.5,0.5];
c=contourcbar;
c.Color='k';
c.FontSize=24;
c.Label.String=['Spice Anomaly'];
ax.Colormap=(customcolormap_preset('pasteljet'));
[conDa,conD]=contour(Months,Depths,CDensInt,[26,25.5,25]);
[conD2a,conD2]=contour(Months,Depths,CDensInt,[26.7,26.3]);
[conD3a,conD3]=contour(Months,Depths,CDensInt,[26.5,27]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-0.6 -0.4 -.2 -0.1 0 0.1 0.2 0.4 0.6];
con.LineColor='none';
conD.LineColor='k';
conD.LineWidth=1.5;
conD.LineStyle='--';
conD.ShowText='on';
conD2.LineColor='k';
conD2.LineWidth=2;
conD2.ShowText='on';
conD3.LineColor='k';
conD3.LineWidth=2;
conD3.ShowText='on';
conD3.LineStyle='-.';
con.LineWidth=0.1;
clabel(conDa,conD,'FontSize',18)
clabel(conD2a,conD2,'FontSize',18)
clabel(conD3a,conD3,'FontSize',18)
ax.YLim=[0,800];
ax.Color=[0.7 0.7 0.7];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=24;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.XTick=min(Months)+1:6:max(Months);
ax.XTickLabel={'Jan 2014','July','Jan 2015','July','Jan 2016','July','Jan 2017','July','Jan 2018'};
ax.Box='on';
ax.LineWidth=1.5;
