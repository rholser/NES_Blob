%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 25-May-2022

%%%%%Required scripts, functions, and toolboxes:
    %sw_dens function (in mixing_library)
    %All_CTD_Data_Bin_V2.mat (from Import_CalcAnomalyBinned_MEOP_V2.m)
%%%%%%%

%Generates T/S plots throughout eseal range in 2x2 degree squares. 
%(randomly removes 1/3 of data points for each plot).

%This script does NOT automatically save all figures generated - figure
%export is done manually.

%Creates map of all data points with rectangles to indicate where T/S plots
%are

%% Load CTD Data
clear
load('All_CTD_Data_Bin_V2.mat');
CTD_Data_Bin.JulDate=ceil(CTD_Data_Bin.JulDate);
CTD_Data_Bin(isnan(CTD_Data_Bin.ConsTempBin),:)=[];

%CTD_Data_Bin=CTD_Data_Bin(CTD_Data_Bin.Depth>=300,:);
depths=unique(CTD_Data_Bin.Depth);
CMean=NaN(size(depths,1),2);

%Upper water masses (0-500m) from "Water Types and Water Masses" in
%Encyclopedia of Ocean Sciences
PSUW=[3 3 15 15 3; 32.6 33.6 33.6 32.6 32.6];
ENPCW=[12 12 20 20 12; 34.2 35.0 35.0 34.2 34.2];
ENPTW=[11 11 20 20 11; 33.8 34.3 34.3 33.8 33.8];

%Intermediat water masses (500-1000m)
PSIW=[5 5 12 12 5; 33.8 34.3 34.3 33.8 33.8];
CIW=[10 10 12 12 10; 33.9 34.4 34.4 33.9 33.9];

NPDW=[1.5 1.5 3 3 1.5; 34.4 34.6 34.6 34.4 34.4];
%% Generate density contours for T-S plots
theta=CTD_Data_Bin.ConsTempBin;
s=CTD_Data_Bin.SalBin;
smin=min(s)-0.01.*min(s);
smax=max(s)+0.01.*max(s);
thetamin=min(theta)-0.1*max(theta);
thetamax=max(theta)+0.1*max(theta);
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
disp(xdim);disp(ydim);
for j=1:ydim
    for i=1:xdim
        dens(j,i)=sw_dens(si(i),thetai(j),0);
    end
end
dens=dens-1000;

%% Define figure parameters
%Axes
Smax=34.5;
Smin=31.5;
Tmax=20;
Tmin=0;
axXC=[0 0 0];
axYC=[0 0 0];
axFSize=20;
axFSizeM= 1.2;
axTC=[0 0 0];

%Density Contours
DStep=0.5;
Space=300;

%Year color and dot size
C14=[0, 0.4470, 0.7410];
C15=[0.8500, 0.3250, 0.0980];
C16=[0.9290, 0.6940, 0.1250];
C17=[0.4940, 0.1840, 0.5560];
GR=[0 0 0];
Dots=10;

%Water Mass Boxes
CPSUW=[0.2, 0.4, 0.2];
CNPDW=[0.6350, 0.0780, 0.1840];
CPSIW=[0, 0.2, 0.6];
CCIW=[];
Lines=2;

%Climatology line/polygon settings
CClim=[0 0 0];
Cpb=[0.6 0.6 0.6];
Cpbe=[0.5 0.5 0.5];
pbLine=1;
pbAlpha=0.5;
spb=0.2; %boundary shrink factor (0-1)

%% Generate all TS Plots

%Create lat/lon structures to cycle through
lat(1,1:2)=[34 36];
lat(2,1:2)=[36 38];
lat(3,1:2)=[38 40];
lat(4,1:2)=[40 42];
lat(5,1:2)=[42 44];
lat(6,1:2)=[44 46];
lat(7,1:2)=[46 48];
lat(8,1:2)=[48 50];
lat(9,1:2)=[50 52];
lat(10,1:2)=[52 54];
lat(11,1:2)=[54 56];
lat(12,1:2)=[56 58];
lat(13,1:2)=[58 60];

lon=NaN(13,31);
lon(1,1:7)=(228:2:240);
lon(2,1:7)=(226:2:238);
lon(3,1:7)=(224:2:236);
lon(4,1:31)=(176:2:236);
lon(5,1:31)=(176:2:236);
lon(6,1:31)=(176:2:236);
lon(7,1:19)=(200:2:236);
lon(8,1:13)=(210:2:234);
lon(9,1:7)=(218:2:230);
lon(10,1:7)=(216:2:228);
lon(11,1:7)=(214:2:226);
lon(12,1:7)=(212:2:224);
lon(13,1:7)=(210:2:222);
cnt2=1;

%Loop through all 2x2 latxlon boxes within the defined range and plot any 
%data present. This will create some empty plots in some cases as not all
%grid cells contain data.
for i=1:13

    data1=CTD_Data_Bin((CTD_Data_Bin.Lat>lat(i,1) & CTD_Data_Bin.Lat<=lat(i,2)),:);
    cnt=1;

    for x=1:5        
        if ~isnan(lon(i,cnt))
            fig=figure(cnt2);
            fig.Color = [1 1 1];
            t=tiledlayout(2,3);
            t.TileSpacing='compact';

            for j=1:6
                data2=data1((data1.Long360<lon(i,cnt+1) & data1.Long360>=lon(i,cnt)),:);
                CMean=NaN(size(depths,1),2);
                for l=1:size(depths,1)
                    CMean(l,1)=mean(data2.CTemp(data2.Depth==depths(l)),'omitnan');
                    CMean(l,2)=mean(data2.CSal(data2.Depth==depths(l)),'omitnan');
                end
                k=boundary(data2.CSal(~isnan(data2.CSal) & ~isnan(data2.CTemp)),...
                    data2.CTemp(~isnan(data2.CSal) & ~isnan(data2.CTemp)),spb);

                nexttile
                ax=gca;
                hold on
                [c,h]=contour(si,thetai,dens,'k');
                h.LevelStep=DStep;
                clabel(c,h,'LabelSpacing',Space);
                ax.YLim=[Tmin,Tmax];
                ax.XLim=[Smin,Smax];
                ax.XColor=axXC;
                ax.YColor=axYC;
                ax.FontSize=axFSize;
                ax.LabelFontSizeMultiplier = axFSizeM;
                ax.Title.Color=axTC;
                if j==1
                    ax.YLabel.String=('\Theta (^oC)');
                elseif j==4
                    ax.YLabel.String=('\Theta (^oC)');
                end
                if j>3
                    ax.XLabel.String='S_{A} (g kg^{-1})';
                end
                ax.Title.String=[num2str(lat(i,1)),'-',num2str(lat(i,2)),'N x ',...
                    num2str(lon(i,cnt)),'-',num2str(lon(i,cnt+1)),'E'];
                if size(data2,1)>0
                    pb=scatter(data2.CSal,data2.CTemp,Dots*2,GR,'filled');
                    int=randi(size(data2,1),round(size(data2,1)/3),1);
                    data2(int(:,1),:)=[];
                    s1=scatter(data2.SalBin(data2.Year==2014),data2.ConsTempBin(data2.Year==2014),...
                        Dots,C14,'filled');
                    s2=scatter(data2.SalBin(data2.Year==2015),data2.ConsTempBin(data2.Year==2015),...
                        Dots,C15,'filled');
                    s3=scatter(data2.SalBin(data2.Year==2016),data2.ConsTempBin(data2.Year==2016),...
                        Dots,C16,'filled');
                    s4=scatter(data2.SalBin(data2.Year==2017),data2.ConsTempBin(data2.Year==2017),...
                        Dots,C17,'filled');
                    p1=plot(PSUW(2,:),PSUW(1,:));
                    p1.LineWidth=Lines;
                    p1.Color=CPSUW;
                    p4=plot(PSIW(2,:),PSIW(1,:));
                    p4.LineWidth=Lines;
                    p4.Color=CPSIW;
                    p=plot(CMean(:,2),CMean(:,1));
                    p.Color=CClim;
                    p.LineWidth=Lines*2;
                if j==4
                    l=legend([s1 s2 s3 s4 pb p1 p4],'2014','2015','2016','2017',...
                        'Climatology','PSUW','PSIW');
                    l.Location='southwest';
                end
                end
                cnt=cnt+1;
                clear k
            end
        else
            cnt=cnt+6;
        end
        cnt2=cnt2+1;
    end
end

%% Map of data and T/S plot pannels
%Indeces for all regions subset by year, salinity data only
Salt=CTD_Data_Bin(~isnan(CTD_Data_Bin.SalBin),:);
ind1=find(Salt.Year==2014);
ind2=find(Salt.Year==2015);
ind3=find(Salt.Year==2016);
ind4=find(Salt.Year==2017);
n=2; %line width

fig=figure(6);
fig.Color = [1 1 1];

g1=geoscatter(Salt.Lat(ind1),Salt.Long360(ind1),20,'filled');
hold on
g2=geoscatter(Salt.Lat(ind2),Salt.Long360(ind2),20,'filled');
g3=geoscatter(Salt.Lat(ind3),Salt.Long360(ind3),20,'filled');
g4=geoscatter(Salt.Lat(ind4),Salt.Long360(ind4),20,'filled');

p=geoplot([34,34,36,36,34],[240,228,228,240,240]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([36,36,38,38,36],[238,226,226,238,238]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([38,38,40,40,38],[236,224,224,236,236]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([40,40,42,42,40],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([42,42,44,44,42],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([44,44,46,46,44],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([46,46,48,48,46],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;

p=geoplot([48,48,50,50,48],[210,222,222,210,210]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([48,48,50,50,48],[222,234,234,222,222]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([50,50,52,52,50],[230,218,218,230,230]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([52,52,54,54,52],[228,216,216,228,228]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([54,54,56,56,54],[226,214,214,226,226]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([56,56,58,58,56],[224,212,212,224,224]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([58,58,60,60,58],[222,210,210,222,222]); 
p.Color='k';
p.LineWidth=n;

geobasemap landcover
ax=gca;
hold off
ax.FontSize=24;
ax.LabelFontSizeMultiplier = 1.2;
ax.Box='on';
ax.LineWidth=1.5;
l=legend([g1 g2 g3 g4],'2014','2015','2016','2017');
ax.Title.String='Data Distribution';

%% Empty Map with Boxes
fig=figure(6);
fig.Color = [1 1 1];

p=geoplot([34,34,36,36,34],[240,228,228,240,240]); 
hold on
p.Color='k';
p.LineWidth=n;

p=geoplot([36,36,38,38,36],[238,226,226,238,238]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([38,38,40,40,38],[236,224,224,236,236]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([40,40,42,42,40],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([40,40,42,42,40],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([42,42,44,44,42],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([42,42,44,44,42],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([44,44,46,46,44],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[188,200,200,188,188]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[176,188,188,176,176]);
p.Color='k';
p.LineWidth=n;

p=geoplot([46,46,48,48,46],[236,224,224,236,236]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[212,224,224,212,212]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[200,212,212,200,200]);
p.Color='k';
p.LineWidth=n;

p=geoplot([48,48,50,50,48],[210,222,222,210,210]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([48,48,50,50,48],[222,234,234,222,222]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([50,50,52,52,50],[230,218,218,230,230]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([52,52,54,54,52],[228,216,216,228,228]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([54,54,56,56,54],[226,214,214,226,226]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([56,56,58,58,56],[224,212,212,224,224]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([58,58,60,60,58],[222,210,210,222,222]); 
p.Color='k';
p.LineWidth=n;

geobasemap grayterrain
ax=gca;
hold off
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.Box='on';
ax.LineWidth=1.5;
l=legend([g1 g2 g3 g4],'2014','2015','2016','2017');
ax.Title.String='2014 - 2017 Data Distribution';