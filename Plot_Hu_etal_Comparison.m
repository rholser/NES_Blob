%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 25-May-2022

%Requires mat file produced by Import_CalcAnomalyBinned_MEOP_V2.m

%Create figures comparable to Hu etal 2017 Fig. 2 (Temp anomaly Depth~Time)
%using seal-collected CTD data for 40-50Nx130-150W.


%% Load CTD Data
clear
load('All_CTD_Data_Bin_V2.mat');
CTD_Data_Bin(isnan(CTD_Data_Bin.ConsTempBin),:)=[];

%% Step 1: subsample data to "Blob" region as defined in Hu etal
CTD_Data_Bin.JulDate=ceil(CTD_Data_Bin.JulDate);

data=CTD_Data_Bin((CTD_Data_Bin.Lat<=50 & CTD_Data_Bin.Lat>=40),:);
data=data((data.Long>=-150 & data.Long<=-130),:);

Dates=(datenum('01/01/2014'):datenum('12/31/2017'))';
Depths=(5:5:1000)';
Months=(1:48)';

%% Step 2: Calculates monthly mean T anomalies at each 5m of depth 
%from 2014-2017, accounting for propogation of error (if data <2sigma,
%assume AnomT is zero)
AnomT=NaN(size(Depths,1),48);
SigmaAnomT=NaN(size(Depths,1),48);
AnomTSigma=NaN(size(Depths,1),48);
TSigma=NaN(size(Depths,1),48);

count=1;
n=2; %threshold number of st deviations to decide to keep data
for i=2014:2017
    casts1=data(data.Year==i,:); 
    for j=1:12
        casts2=casts1(casts1.Month==j,:);
        for k=1:size(Depths,1)-1
            nT=size(casts2.AnomConsT(casts2.Depth==Depths(k)),1);
            AnomT(k,count)=mean(casts2.AnomConsT(casts2.Depth==Depths(k)),'omitnan');
            AnomTSigma(k,count)=mean(casts2.AnomConsTSigmas(casts2.Depth==Depths(k)),'omitnan');
            TSigma(k,count)=mean(casts2.SigmaAnomT(casts2.Depth==Depths(k)),'omitnan');
            SigmaAnomT(k,count) = sqrt((1/nT^2)*sum((casts2.SigmaAnomT(casts2.Depth==Depths(k)).^2),'omitnan'));
            if abs(AnomT(k,count))<(n*SigmaAnomT(k,count))
                AnomT(k,count)=0;
            end
        end
        count=count+1;
    end  
    clear casts1 casts2
end

%% Step 3: Calculates depth interpolation of means from Step 2 above.
AnomT2=NaN(size(Depths,1),48);
AnomTSigma2=NaN(size(Depths,1),48);
TSigma2=NaN(size(Depths,1),48);

for i=1:48
    ind1=find(~isnan(AnomT(:,i)));
    if size(ind1,1)>0
        d1=Depths(ind1);
        data1=AnomT(ind1,i);
        AnomT2(:,i)=interp1(d1,data1,Depths);
        data2=AnomTSigma(ind1,i);
        AnomTSigma2(:,i)=interp1(d1,data2,Depths);
        data3=TSigma(ind1,i);
        TSigma2(:,i)=interp1(d1,data3,Depths);
        
    end
    clear ind1 d1 data1 data2
end

%% Step 4: Plot interpolated data contours
%Temperature Anomalies
fig=figure(1);
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Months,Depths,AnomT2,[-2.5 -2 -1.5 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1 1.5 2 2.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-2 -1 -0.5 -.2 0 0.2 0.5 1 2];
con.LineColor='k';
con.LineWidth=0.5;
clabel(con1,con,'FontSize',18)
c=contourcbar;
c.Color='k';
c.FontSize=24;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.CLim=[-2.5,2.5];
ax.YLim=[0,800];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=24;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.XTick=min(Months)+1:6:max(Months);
ax.XTickLabel={'Jan','July','Jan','July','Jan','July','Jan','July','Jan'};
ax.Box='on';
ax.LineWidth=1.5;

fig=figure(2);
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Months,Depths,AnomTSigma2,...
    [-2.5 -2 -1.5 -0.5 0 0.5 1 1.5 2 3 4 5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-2 -1 -0.5 -.01 0 0.1 0.5 1 2 3 4 5];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',18)
c=colorbar;
c.Color='k';
c.FontSize=24;
c.Label.String=['Temperature Anomaly/Sigma'];
ax.Colormap=([0 0 0;redblue]);
ax.CLim=[-5,5];
ax.YLim=[0,800];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=24;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.XTick=min(Months)+1:6:max(Months);
ax.XTickLabel={'Jan','July','Jan','July','Jan','July','Jan','July','Jan'};
ax.Box='on';
ax.LineWidth=1.5;

fig=figure(3);
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Months,Depths,TSigma2,...
    [0:0.1:1.6]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[0:0.1:1.6];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',16,'Color','w')
c=contourcbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Sigma T'];
ax.Colormap=parula;
ax.CLim=[0,1.5];
ax.YLim=[0,800];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=22;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.XTick=min(Months)+1:6:max(Months);
ax.XTickLabel={'Jan','July','Jan','July','Jan','July','Jan','July','Jan'};
ax.Box='on';
ax.LineWidth=1.5;
