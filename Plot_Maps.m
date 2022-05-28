%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 25-May-2022

%%%%%Required scripts, functions, and toolboxes:
    %All_CTD_Data_Bin_V2.mat (from Import_CalcAnomalyBinned_MEOP_V2.m)
%%%%%%%
%Generates geoaxis maps of CTD data locations with ocean/land background
%and overlays regions/subdivisions of the data

%Script does not save figures, export must be done manually.


%% Load Data
clear

%CTD data
load('All_CTD_Data_Bin_V2.mat');
CTD_Data_Bin.JulDate=ceil(CTD_Data_Bin.JulDate);

%Indeces for all regions subset by year
ind1=find(CTD_Data_Bin.Year==2014);
ind2=find(CTD_Data_Bin.Year==2015);
ind3=find(CTD_Data_Bin.Year==2016);
ind4=find(CTD_Data_Bin.Year==2017);

n=2; %line width
%% Map of Blob rectangles and year
fig=figure(1);
fig.Color = [1 1 1];

g1=geoscatter(CTD_Data_Bin.Lat(ind1),CTD_Data_Bin.Long360(ind1),10,'filled');
hold on
g2=geoscatter(CTD_Data_Bin.Lat(ind2),CTD_Data_Bin.Long360(ind2),10,'filled');
g3=geoscatter(CTD_Data_Bin.Lat(ind3),CTD_Data_Bin.Long360(ind3),10,'filled');
g4=geoscatter(CTD_Data_Bin.Lat(ind4),CTD_Data_Bin.Long360(ind4),10,'filled');

p=geoplot([40,40,50,50,40],[230,210,210,230,230]); 
p.Color='k';
p.LineWidth=n;

p=geoplot([40,40,45,45,40],[236,180,180,236,236]); 
p.Color='k';
p.LineWidth=n;
p.LineStyle='--';

geobasemap landcover
ax=gca;
hold off
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.Box='on';
ax.LineWidth=1.5;
ax.Title.String='Data Distribution';
l=legend([g1 g2 g3 g4],'2014','2015','2016','2017');



%% Map of T/S location squares and year
fig=figure(1);
fig.Color = [1 1 1];

g1=geoscatter(CTD_Data_Bin.Lat(ind1),CTD_Data_Bin.Long360(ind1),10,'filled');
hold on
g2=geoscatter(CTD_Data_Bin.Lat(ind2),CTD_Data_Bin.Long360(ind2),10,'filled');
g3=geoscatter(CTD_Data_Bin.Lat(ind3),CTD_Data_Bin.Long360(ind3),10,'filled');
g4=geoscatter(CTD_Data_Bin.Lat(ind4),CTD_Data_Bin.Long360(ind4),10,'filled');

p=geoplot([44,44,46,46,44],[226,224,224,226,226]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[228,226,226,228,228]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[230,228,228,230,230]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[232,230,230,232,232]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[234,232,232,234,234]);
p.Color='k';
p.LineWidth=n;
p=geoplot([44,44,46,46,44],[236,234,234,236,236]);
p.Color='k';
p.LineWidth=n;

p=geoplot([46,46,48,48,46],[212,214,214,212,212]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[214,216,216,214,214]); 
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[216,218,218,216,216]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[218,220,220,218,218]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[220,222,222,220,220]);
p.Color='k';
p.LineWidth=n;
p=geoplot([46,46,48,48,46],[222,224,224,222,222]);
p.Color='k';
p.LineWidth=n;

geobasemap landcover
ax=gca;
hold off
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.Box='on';
ax.LineWidth=1.5;
l=legend([g1 g2 g3 g4],'2014','2015','2016','2017');