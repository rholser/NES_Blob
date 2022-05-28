%Created by: Rachel Holser (rholser@ucsc.edu)
%Last Updated: 25-May-2021

%Requires mat file produced by Import_CalcAnomalyBinned_MEOP_V2.m

%Creates Longitude sections of temperature, salinity, and density anomalies of
%seasonal (3-month) means from 2014-2017 (each figure is a 4x4 grid)

% Step 1 - Calculate mean Anomaly values at each depth bin(5m)/longitude(1deg)
% Step 2 - Calculates sigmaTa for each depth layer when possible, using sigmaTtag 
%           (0.02 - Roquet et al. 2011) and sigmaTc. 
% Step 3 - Excludes AnomT within 2 SigmaTa.
% Step 4 - Linear interpolation of AnomT values to fill depth gaps

%Generates tiled figures of calculated data

%% Load only CTD Data
clear
load('All_CTD_Data_Bin_V2.mat');

%%Depth~Longitude of values from Y1 to Y2 degrees N, 3-month means

%Subset data to Latitudes Y1 and Y2
Y1=40;
Y2=44;
data=CTD_Data_Bin((CTD_Data_Bin.Lat<Y2 & CTD_Data_Bin.Lat>=Y1),:);
data.Long=round(data.Long);

%Subsest data to 3-month intervals
data1=data((data.JulDate>=datenum('12-01-2013','mm-dd-yyyy') & ...
         data.JulDate<=datenum('02-28-2014','mm-dd-yyyy')),:);
data2=data((data.JulDate>=datenum('03-01-2014','mm-dd-yyyy') & ...
         data.JulDate<=datenum('05-31-2014','mm-dd-yyyy')),:);
data3=data((data.JulDate>=datenum('06-01-2014','mm-dd-yyyy') & ...
         data.JulDate<=datenum('08-31-2014','mm-dd-yyyy')),:);
data4=data((data.JulDate>=datenum('09-01-2014','mm-dd-yyyy') & ...
         data.JulDate<=datenum('11-30-2014','mm-dd-yyyy')),:);
data5=data((data.JulDate>=datenum('12-01-2014','mm-dd-yyyy') & ...
         data.JulDate<=datenum('02-28-2015','mm-dd-yyyy')),:);
data6=data((data.JulDate>=datenum('03-01-2015','mm-dd-yyyy') & ...
         data.JulDate<=datenum('05-31-2015','mm-dd-yyyy')),:);
data7=data((data.JulDate>=datenum('06-01-2015','mm-dd-yyyy') & ...
         data.JulDate<=datenum('08-31-2015','mm-dd-yyyy')),:);
data8=data((data.JulDate>=datenum('09-01-2015','mm-dd-yyyy') & ...
         data.JulDate<=datenum('11-30-2015','mm-dd-yyyy')),:);   
data9=data((data.JulDate>=datenum('12-01-2015','mm-dd-yyyy') & ...
         data.JulDate<=datenum('02-28-2016','mm-dd-yyyy')),:);
data10=data((data.JulDate>=datenum('03-01-2016','mm-dd-yyyy') & ...
         data.JulDate<=datenum('05-31-2016','mm-dd-yyyy')),:);
data11=data((data.JulDate>=datenum('06-01-2016','mm-dd-yyyy') & ...
         data.JulDate<=datenum('08-31-2016','mm-dd-yyyy')),:);
data12=data((data.JulDate>=datenum('09-01-2016','mm-dd-yyyy') & ...
         data.JulDate<=datenum('11-30-2016','mm-dd-yyyy')),:);
data13=data((data.JulDate>=datenum('12-01-2016','mm-dd-yyyy') & ...
         data.JulDate<=datenum('02-28-2017','mm-dd-yyyy')),:);
data14=data((data.JulDate>=datenum('03-01-2017','mm-dd-yyyy') & ...
         data.JulDate<=datenum('05-31-2017','mm-dd-yyyy')),:);
data15=data((data.JulDate>=datenum('06-01-2017','mm-dd-yyyy') & ...
         data.JulDate<=datenum('08-31-2017','mm-dd-yyyy')),:);
data16=data((data.JulDate>=datenum('09-01-2017','mm-dd-yyyy') & ...
         data.JulDate<=datenum('11-30-2017','mm-dd-yyyy')),:);

%% Calculate mean temperature anomaly values for Longitude and 5m depth bins
Longitudes=(-180:-124)'; %Longitude range used
Depths=(5:5:1000)'; %resolution for interpolation

%Preallocate structures for calculated values - temperature
AnomT1=NaN(size(Depths,1),size(Longitudes,1));
AnomT2=NaN(size(Depths,1),size(Longitudes,1));
AnomT3=NaN(size(Depths,1),size(Longitudes,1));
AnomT4=NaN(size(Depths,1),size(Longitudes,1));
AnomT5=NaN(size(Depths,1),size(Longitudes,1));
AnomT6=NaN(size(Depths,1),size(Longitudes,1));
AnomT7=NaN(size(Depths,1),size(Longitudes,1));
AnomT8=NaN(size(Depths,1),size(Longitudes,1));
AnomT9=NaN(size(Depths,1),size(Longitudes,1));
AnomT10=NaN(size(Depths,1),size(Longitudes,1));
AnomT11=NaN(size(Depths,1),size(Longitudes,1));
AnomT12=NaN(size(Depths,1),size(Longitudes,1));
AnomT13=NaN(size(Depths,1),size(Longitudes,1));
AnomT14=NaN(size(Depths,1),size(Longitudes,1));
AnomT15=NaN(size(Depths,1),size(Longitudes,1));
AnomT16=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT1=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT2=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT3=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT4=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT5=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT6=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT7=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT8=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT9=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT10=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT11=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT12=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT13=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT14=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT15=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomT16=NaN(size(Depths,1),size(Longitudes,1));
%Preallocate structures for calculated values - salinity
AnomS1=NaN(size(Depths,1),size(Longitudes,1));
AnomS2=NaN(size(Depths,1),size(Longitudes,1));
AnomS3=NaN(size(Depths,1),size(Longitudes,1));
AnomS4=NaN(size(Depths,1),size(Longitudes,1));
AnomS5=NaN(size(Depths,1),size(Longitudes,1));
AnomS6=NaN(size(Depths,1),size(Longitudes,1));
AnomS7=NaN(size(Depths,1),size(Longitudes,1));
AnomS8=NaN(size(Depths,1),size(Longitudes,1));
AnomS9=NaN(size(Depths,1),size(Longitudes,1));
AnomS10=NaN(size(Depths,1),size(Longitudes,1));
AnomS11=NaN(size(Depths,1),size(Longitudes,1));
AnomS12=NaN(size(Depths,1),size(Longitudes,1));
AnomS13=NaN(size(Depths,1),size(Longitudes,1));
AnomS14=NaN(size(Depths,1),size(Longitudes,1));
AnomS15=NaN(size(Depths,1),size(Longitudes,1));
AnomS16=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS1=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS2=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS3=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS4=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS5=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS6=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS7=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS8=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS9=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS10=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS11=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS12=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS13=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS14=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS15=NaN(size(Depths,1),size(Longitudes,1));
SigmaAnomS16=NaN(size(Depths,1),size(Longitudes,1));
%Preallocate structures for calculated values - density
AnomD1=NaN(size(Depths,1),size(Longitudes,1));
AnomD2=NaN(size(Depths,1),size(Longitudes,1));
AnomD3=NaN(size(Depths,1),size(Longitudes,1));
AnomD4=NaN(size(Depths,1),size(Longitudes,1));
AnomD5=NaN(size(Depths,1),size(Longitudes,1));
AnomD6=NaN(size(Depths,1),size(Longitudes,1));
AnomD7=NaN(size(Depths,1),size(Longitudes,1));
AnomD8=NaN(size(Depths,1),size(Longitudes,1));
AnomD9=NaN(size(Depths,1),size(Longitudes,1));
AnomD10=NaN(size(Depths,1),size(Longitudes,1));
AnomD11=NaN(size(Depths,1),size(Longitudes,1));
AnomD12=NaN(size(Depths,1),size(Longitudes,1));
AnomD13=NaN(size(Depths,1),size(Longitudes,1));
AnomD14=NaN(size(Depths,1),size(Longitudes,1));
AnomD15=NaN(size(Depths,1),size(Longitudes,1));
AnomD16=NaN(size(Depths,1),size(Longitudes,1));
D1=NaN(size(Depths,1),size(Longitudes,1));
D2=NaN(size(Depths,1),size(Longitudes,1));
D3=NaN(size(Depths,1),size(Longitudes,1));
D4=NaN(size(Depths,1),size(Longitudes,1));
D5=NaN(size(Depths,1),size(Longitudes,1));
D6=NaN(size(Depths,1),size(Longitudes,1));
D7=NaN(size(Depths,1),size(Longitudes,1));
D8=NaN(size(Depths,1),size(Longitudes,1));
D9=NaN(size(Depths,1),size(Longitudes,1));
D10=NaN(size(Depths,1),size(Longitudes,1));
D11=NaN(size(Depths,1),size(Longitudes,1));
D12=NaN(size(Depths,1),size(Longitudes,1));
D13=NaN(size(Depths,1),size(Longitudes,1));
D14=NaN(size(Depths,1),size(Longitudes,1));
D15=NaN(size(Depths,1),size(Longitudes,1));
D16=NaN(size(Depths,1),size(Longitudes,1));

n=2; %threshold number of st deviations within which to keep data

%Sesonal calculation of mean and sigma anomaly
for i=1:size(Longitudes,1)-1
    %for each 3-month interval, find data at Longitude(i)
    casts1=data1(data1.Long==Longitudes(i),:);
    casts2=data2(data2.Long==Longitudes(i),:);
    casts3=data3(data3.Long==Longitudes(i),:);
    casts4=data4(data4.Long==Longitudes(i),:);
    casts5=data5(data5.Long==Longitudes(i),:);
    casts6=data6(data6.Long==Longitudes(i),:);
    casts7=data7(data7.Long==Longitudes(i),:);
    casts8=data8(data8.Long==Longitudes(i),:);
    casts9=data9(data9.Long==Longitudes(i),:);
    casts10=data10(data10.Long==Longitudes(i),:);
    casts11=data11(data11.Long==Longitudes(i),:);
    casts12=data12(data12.Long==Longitudes(i),:);
    casts13=data13(data13.Long==Longitudes(i),:);
    casts14=data14(data14.Long==Longitudes(i),:);
    casts15=data15(data15.Long==Longitudes(i),:);
    casts16=data16(data16.Long==Longitudes(i),:);

    %for each Depth(k) caculate mean and sigma of anomalies in T, S, and D.
    %Exclude data greater than 2 sigmas from the mean
    for k=1:size(Depths,1)-1
        temp=casts1(casts1.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT1(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT1(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT1(k,i))<(n*SigmaAnomT1(k,i))
                AnomT1(k,i)=0;
            end
            SigmaAnomS1(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS1(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS1(k,i))<(n*SigmaAnomS1(k,i))
                AnomS1(k,i)=0;
            end
            AnomD1(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D1(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts2(casts2.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT2(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT2(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT2(k,i))<(n*SigmaAnomT2(k,i))
                AnomT2(k,i)=0;
            end
            SigmaAnomS2(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS2(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS2(k,i))<(n*SigmaAnomS2(k,i))
                AnomS2(k,i)=0;
            end
            AnomD2(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D2(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts3(casts3.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT3(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT3(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT3(k,i))<(n*SigmaAnomT3(k,i))
                AnomT3(k,i)=0;
            end
            SigmaAnomS3(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS3(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS3(k,i))<(n*SigmaAnomS3(k,i))
                AnomS3(k,i)=0;
            end
            AnomD3(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D3(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts4(casts4.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT4(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT4(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT4(k,i))<(n*SigmaAnomT4(k,i))
                AnomT4(k,i)=0;
            end
            SigmaAnomS4(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS4(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS4(k,i))<(n*SigmaAnomS4(k,i))
                AnomS4(k,i)=0;
            end
            AnomD4(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D4(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts5(casts5.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT5(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT5(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT5(k,i))<(n*SigmaAnomT5(k,i))
                AnomT5(k,i)=0;
            end
            SigmaAnomS5(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS5(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS5(k,i))<(n*SigmaAnomS5(k,i))
                AnomS5(k,i)=0;
            end
            AnomD5(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D5(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts6(casts6.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT6(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT6(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT6(k,i))<(n*SigmaAnomT6(k,i))
                AnomT6(k,i)=0;
            end
            SigmaAnomS6(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS6(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS6(k,i))<(n*SigmaAnomS6(k,i))
                AnomS6(k,i)=0;
            end
            AnomD6(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D6(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts7(casts7.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT7(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT7(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT7(k,i))<(n*SigmaAnomT7(k,i))
                AnomT7(k,i)=0;
            end
            SigmaAnomS7(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS7(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS7(k,i))<(n*SigmaAnomS7(k,i))
                AnomS7(k,i)=0;
            end
            AnomD7(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D7(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts8(casts8.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT8(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT8(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT8(k,i))<(n*SigmaAnomT8(k,i))
                AnomT8(k,i)=0;
            end
            SigmaAnomS8(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS8(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS8(k,i))<(n*SigmaAnomS8(k,i))
                AnomS8(k,i)=0;
            end
            AnomD8(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D8(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts9(casts9.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT9(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT9(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT9(k,i))<(n*SigmaAnomT9(k,i))
                AnomT9(k,i)=0;
            end
            SigmaAnomS9(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS9(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS9(k,i))<(n*SigmaAnomS9(k,i))
                AnomS9(k,i)=0;
            end
            AnomD9(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D9(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts10(casts10.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT10(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT10(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT10(k,i))<(n*SigmaAnomT10(k,i))
                AnomT10(k,i)=0;
            end
            SigmaAnomS10(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS10(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS10(k,i))<(n*SigmaAnomS10(k,i))
                AnomS10(k,i)=0;
            end
            AnomD10(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D10(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts11(casts11.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT11(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT11(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT11(k,i))<(n*SigmaAnomT11(k,i))
                AnomT11(k,i)=0;
            end
            SigmaAnomS11(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS11(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS11(k,i))<(n*SigmaAnomS11(k,i))
                AnomS11(k,i)=0;
            end
            AnomD11(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D11(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts12(casts12.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT12(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT12(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT12(k,i))<(n*SigmaAnomT12(k,i))
                AnomT12(k,i)=0;
            end
            SigmaAnomS12(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS12(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS12(k,i))<(n*SigmaAnomS12(k,i))
                AnomS12(k,i)=0;
            end
            AnomD12(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D12(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts13(casts13.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT13(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT13(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT13(k,i))<(n*SigmaAnomT13(k,i))
                AnomT13(k,i)=0;
            end
            SigmaAnomS13(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS13(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS13(k,i))<(n*SigmaAnomS13(k,i))
                AnomS13(k,i)=0;
            end
            AnomD13(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D13(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts14(casts14.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT14(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT14(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT14(k,i))<(n*SigmaAnomT14(k,i))
                AnomT14(k,i)=0;
            end
            SigmaAnomS14(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS14(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS14(k,i))<(n*SigmaAnomS14(k,i))
                AnomS14(k,i)=0;
            end
            AnomD14(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D14(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts15(casts15.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT1(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT15(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT15(k,i))<(n*SigmaAnomT15(k,i))
                AnomT15(k,i)=0;
            end
            SigmaAnomS15(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS15(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS15(k,i))<(n*SigmaAnomS15(k,i))
                AnomS15(k,i)=0;
            end
            AnomD15(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D15(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
        temp=casts16(casts16.Depth==Depths(k),:);
        if size(temp,1)>0
            SigmaAnomT16(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomT.^2,'omitnan'));
            AnomT16(k,i)=mean(temp.AnomConsT(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomT16(k,i))<(n*SigmaAnomT16(k,i))
                AnomT16(k,i)=0;
            end
            SigmaAnomS16(k,i) = sqrt((1/(size(temp,1)^2))*sum(temp.SigmaAnomS.^2,'omitnan'));
            AnomS16(k,i)=mean(temp.AnomAbsSal(temp.Depth==Depths(k)),'omitnan');
            if abs(AnomS16(k,i))<(n*SigmaAnomS16(k,i))
                AnomS16(k,i)=0;
            end
            AnomD16(k,i)=mean(temp.AnomD(temp.Depth==Depths(k)),'omitnan');
            D16(k,i)=mean(temp.Density(temp.Depth==Depths(k)),'omitnan');
        end
    end
        clear casts1 casts2 casts3 casts4 casts5 casts6 casts7 casts8 casts9...
        casts10 casts11 casts12 casts13 casts14 casts15 casts16
end

%Linear interpolation of all calculated values to standard 5m depth
%resolution
for i=1:size(Longitudes,1)
    ind1=find(~isnan(AnomT1(:,i)));
    ind2=find(~isnan(AnomT2(:,i)));
    ind3=find(~isnan(AnomT3(:,i)));
    ind4=find(~isnan(AnomT4(:,i)));
    ind5=find(~isnan(AnomT5(:,i)));
    ind6=find(~isnan(AnomT6(:,i)));
    ind7=find(~isnan(AnomT7(:,i)));
    ind8=find(~isnan(AnomT8(:,i)));
    ind9=find(~isnan(AnomT9(:,i)));
    ind10=find(~isnan(AnomT10(:,i)));
    ind11=find(~isnan(AnomT11(:,i)));
    ind12=find(~isnan(AnomT12(:,i)));
    ind13=find(~isnan(AnomT13(:,i)));
    ind14=find(~isnan(AnomT14(:,i)));
    ind15=find(~isnan(AnomT15(:,i)));
    ind16=find(~isnan(AnomT16(:,i)));
    ind1S=find(~isnan(AnomS1(:,i)));
    ind2S=find(~isnan(AnomS2(:,i)));
    ind3S=find(~isnan(AnomS3(:,i)));
    ind4S=find(~isnan(AnomS4(:,i)));
    ind5S=find(~isnan(AnomS5(:,i)));
    ind6S=find(~isnan(AnomS6(:,i)));
    ind7S=find(~isnan(AnomS7(:,i)));
    ind8S=find(~isnan(AnomS8(:,i)));
    ind9S=find(~isnan(AnomS9(:,i)));
    ind10S=find(~isnan(AnomS10(:,i)));
    ind11S=find(~isnan(AnomS11(:,i)));
    ind12S=find(~isnan(AnomS12(:,i)));
    ind13S=find(~isnan(AnomS13(:,i)));
    ind14S=find(~isnan(AnomS14(:,i)));
    ind15S=find(~isnan(AnomS15(:,i)));
    ind16S=find(~isnan(AnomS16(:,i)));
    ind1D=find(~isnan(AnomD1(:,i)));
    ind2D=find(~isnan(AnomD2(:,i)));
    ind3D=find(~isnan(AnomD3(:,i)));
    ind4D=find(~isnan(AnomD4(:,i)));
    ind5D=find(~isnan(AnomD5(:,i)));
    ind6D=find(~isnan(AnomD6(:,i)));
    ind7D=find(~isnan(AnomD7(:,i)));
    ind8D=find(~isnan(AnomD8(:,i)));
    ind9D=find(~isnan(AnomD9(:,i)));
    ind10D=find(~isnan(AnomD10(:,i)));
    ind11D=find(~isnan(AnomD11(:,i)));
    ind12D=find(~isnan(AnomD12(:,i)));
    ind13D=find(~isnan(AnomD13(:,i)));
    ind14D=find(~isnan(AnomD14(:,i)));
    ind15D=find(~isnan(AnomD15(:,i)));
    ind16D=find(~isnan(AnomD16(:,i)));
    indD1=find(~isnan(D1(:,i)));
    indD2=find(~isnan(D2(:,i)));
    indD3=find(~isnan(D3(:,i)));
    indD4=find(~isnan(D4(:,i)));
    indD5=find(~isnan(D5(:,i)));
    indD6=find(~isnan(D6(:,i)));
    indD7=find(~isnan(D7(:,i)));
    indD8=find(~isnan(D8(:,i)));
    indD9=find(~isnan(D9(:,i)));
    indD10=find(~isnan(D10(:,i)));
    indD11=find(~isnan(D11(:,i)));
    indD12=find(~isnan(D12(:,i)));
    indD13=find(~isnan(D13(:,i)));
    indD14=find(~isnan(D14(:,i)));
    indD15=find(~isnan(D15(:,i)));
    indD16=find(~isnan(D16(:,i)));
    
    if size(ind1,1)>1
        d1=Depths(ind1);
        data1=AnomT1(ind1,i);
        AnomT1(:,i)=interp1(d1,data1,Depths);
    end
    if size(ind1S,1)>1
        d1S=Depths(ind1S);
        data1S=AnomS1(ind1S,i);
        AnomS1(:,i)=interp1(d1S,data1S,Depths);
    end
    if size(ind1D,1)>1
        d1D=Depths(ind1D);
        data1D=AnomD1(ind1D,i);
        AnomD1(:,i)=interp1(d1D,data1D,Depths);
    end
    if size(indD1,1)>1
        dD1=Depths(indD1);
        dataD1=D1(indD1,i);
        D1(:,i)=interp1(dD1,dataD1,Depths);
    end
    if size(ind2,1)>1
        d2=Depths(ind2);
        data2=AnomT2(ind2,i);
        AnomT2(:,i)=interp1(d2,data2,Depths);
    end
    if size(ind2S,1)>1
        d2S=Depths(ind2S);
        data2S=AnomS2(ind2S,i);
        AnomS2(:,i)=interp1(d2S,data2S,Depths);
    end
    if size(ind2D,1)>1
        d2D=Depths(ind2D);
        data2D=AnomD2(ind2D,i);
        AnomD2(:,i)=interp1(d2D,data2D,Depths);
    end
    if size(indD2,1)>1
        dD2=Depths(indD2);
        dataD2=D2(indD2,i);
        D2(:,i)=interp1(dD2,dataD2,Depths);
    end
    
    if size(ind3,1)>1
        d3=Depths(ind3);
        data3=AnomT3(ind3,i);
        AnomT3(:,i)=interp1(d3,data3,Depths);
    end
    if size(ind3S,1)>1
        d3S=Depths(ind3S);
        data3S=AnomS3(ind3S,i);
        AnomS3(:,i)=interp1(d3S,data3S,Depths);
    end
    if size(ind3D,1)>1
        d3D=Depths(ind3D);
        data3D=AnomD3(ind3D,i);
        AnomD3(:,i)=interp1(d3D,data3D,Depths);
    end
    if size(indD3,1)>1
        dD3=Depths(indD3);
        dataD3=D3(indD3,i);
        D3(:,i)=interp1(dD3,dataD3,Depths);
    end
    if size(ind4,1)>1
        d4=Depths(ind4);
        data4=AnomT4(ind4,i);
        AnomT4(:,i)=interp1(d4,data4,Depths);
    end
    if size(ind4S,1)>1
        d4S=Depths(ind4S);
        data4S=AnomS4(ind4S,i);
        AnomS4(:,i)=interp1(d4S,data4S,Depths);
    end
    if size(ind4D,1)>1
        d4D=Depths(ind4D);
        data4D=AnomS4(ind4D,i);
        AnomD4(:,i)=interp1(d4D,data4D,Depths);
    end
    if size(indD4,1)>1
        dD4=Depths(indD4);
        dataD4=D4(indD4,i);
        D4(:,i)=interp1(dD4,dataD4,Depths);
    end
    if size(ind5,1)>1
        d5=Depths(ind5);
        data5=AnomT5(ind5,i);
        AnomT5(:,i)=interp1(d5,data5,Depths);
    end
    if size(ind5S,1)>1
        d5S=Depths(ind5S);
        data5S=AnomS5(ind5S,i);
        AnomS5(:,i)=interp1(d5S,data5S,Depths);
    end
    if size(ind5D,1)>1
        d5D=Depths(ind5D);
        data5D=AnomD5(ind5D,i);
        AnomD5(:,i)=interp1(d5D,data5D,Depths);
    end
    if size(indD5,1)>1
        dD5=Depths(indD5);
        dataD5=D5(indD5,i);
        D5(:,i)=interp1(dD5,dataD5,Depths);
    end
    if size(ind6,1)>1
        d6=Depths(ind6);
        data6=AnomT6(ind6,i);
        AnomT6(:,i)=interp1(d6,data6,Depths);
    end
    if size(ind6S,1)>1
        d6S=Depths(ind6S);
        data6S=AnomS6(ind6S,i);
        AnomS6(:,i)=interp1(d6S,data6S,Depths);
    end
    if size(ind6D,1)>1
        d6D=Depths(ind6D);
        data6D=AnomD6(ind6D,i);
        AnomD6(:,i)=interp1(d6D,data6D,Depths);
    end
    if size(indD6,1)>1
        dD6=Depths(indD6);
        dataD6=D6(indD6,i);
        D6(:,i)=interp1(dD6,dataD6,Depths);
    end
    if size(ind7,1)>1
        d7=Depths(ind7);
        data7=AnomT7(ind7,i);
        AnomT7(:,i)=interp1(d7,data7,Depths);
    end
    if size(ind7S,1)>1
        d7S=Depths(ind7S);
        data7S=AnomS7(ind7S,i);
        AnomS7(:,i)=interp1(d7S,data7S,Depths);
    end
    if size(ind7D,1)>1
        d7D=Depths(ind7D);
        data7D=AnomD7(ind7D,i);
        AnomD7(:,i)=interp1(d7D,data7D,Depths);
    end
    if size(indD7,1)>1
        dD7=Depths(indD7);
        dataD7=D7(indD7,i);
        D7(:,i)=interp1(dD7,dataD7,Depths);
    end
    if size(ind8,1)>1
        d8=Depths(ind8);
        data8=AnomT8(ind8,i);
        AnomT8(:,i)=interp1(d8,data8,Depths);
    end
    if size(ind8S,1)>1
        d8S=Depths(ind8S);
        data8S=AnomS8(ind8S,i);
        AnomS8(:,i)=interp1(d8S,data8S,Depths);
    end
    if size(ind8D,1)>1
        d8D=Depths(ind8D);
        data8D=AnomD8(ind8D,i);
        AnomD8(:,i)=interp1(d8D,data8D,Depths);
    end
    if size(indD8,1)>1
        dD8=Depths(indD8);
        dataD8=D8(indD8,i);
        D8(:,i)=interp1(dD8,dataD8,Depths);
    end
    if size(ind9,1)>1
        d9=Depths(ind9);
        data9=AnomT9(ind9,i);
        AnomT9(:,i)=interp1(d9,data9,Depths);
    end
    if size(ind9S,1)>1
        d9S=Depths(ind9S);
        data9S=AnomS9(ind9S,i);
        AnomS9(:,i)=interp1(d9S,data9S,Depths);
    end
    if size(ind9D,1)>1
        d9D=Depths(ind9D);
        data9D=AnomD9(ind9D,i);
        AnomD9(:,i)=interp1(d9D,data9D,Depths);
    end
    if size(indD9,1)>1
        dD9=Depths(indD9);
        dataD9=D9(indD9,i);
        D9(:,i)=interp1(dD9,dataD9,Depths);
    end
    if size(ind10,1)>1
        d10=Depths(ind10);
        data10=AnomT10(ind10,i);
        AnomT10(:,i)=interp1(d10,data10,Depths);
    end
    if size(ind10S,1)>1
        d10S=Depths(ind10S);
        data10S=AnomS10(ind10S,i);
        AnomS10(:,i)=interp1(d10S,data10S,Depths);
    end
    if size(ind10D,1)>1
        d10D=Depths(ind10D);
        data10D=AnomD10(ind10D,i);
        AnomD10(:,i)=interp1(d10D,data10D,Depths);
    end
    if size(indD10,1)>1
        dD10=Depths(indD10);
        dataD10=D10(indD10,i);
        D10(:,i)=interp1(dD10,dataD10,Depths);
    end
    if size(ind11,1)>1
        d11=Depths(ind11);
        data11=AnomT11(ind11,i);
        AnomT11(:,i)=interp1(d11,data11,Depths);
    end
    if size(ind11S,1)>1
        d11S=Depths(ind11S);
        data11S=AnomS11(ind11S,i);
        AnomS11(:,i)=interp1(d11S,data11S,Depths);
    end
    if size(ind11D,1)>1
        d11D=Depths(ind11D);
        data11D=AnomD11(ind11D,i);
        AnomD11(:,i)=interp1(d11D,data11D,Depths);
    end
    if size(indD11,1)>1
        dD11=Depths(indD11);
        dataD11=D11(indD11,i);
        D11(:,i)=interp1(dD11,dataD11,Depths);
    end
    if size(ind12,1)>1
        d12=Depths(ind12);
        data12=AnomT12(ind12,i);
        AnomT12(:,i)=interp1(d12,data12,Depths);
    end
    if size(ind12S,1)>1
        d12S=Depths(ind12S);
        data12S=AnomS12(ind12S,i);
        AnomS12(:,i)=interp1(d12S,data12S,Depths);
    end
    if size(ind12D,1)>1
        d12D=Depths(ind12D);
        data12D=AnomD12(ind12D,i);
        AnomD12(:,i)=interp1(d12D,data12D,Depths);
    end
    if size(indD12,1)>1
        dD12=Depths(indD12);
        dataD12=D12(indD12,i);
        D12(:,i)=interp1(dD12,dataD12,Depths);
    end
    if size(ind13,1)>1
        d13=Depths(ind13);
        data13=AnomT13(ind13,i);
        AnomT13(:,i)=interp1(d13,data13,Depths);
    end
    if size(ind13S,1)>1
        d13S=Depths(ind13S);
        data13S=AnomS13(ind13S,i);
        AnomS13(:,i)=interp1(d13S,data13S,Depths);
    end
    if size(ind13D,1)>1
        d13D=Depths(ind13D);
        data13D=AnomD13(ind13D,i);
        AnomD13(:,i)=interp1(d13D,data13D,Depths);
    end
    if size(indD13,1)>1
        dD13=Depths(indD13);
        dataD13=D13(indD13,i);
        D13(:,i)=interp1(dD13,dataD13,Depths);
    end
    if size(ind14,1)>1
        d14=Depths(ind14);
        data14=AnomT14(ind14,i);
        AnomT14(:,i)=interp1(d14,data14,Depths);
    end
    if size(ind14S,1)>1
        d14S=Depths(ind14S);
        data14S=AnomS14(ind14S,i);
        AnomS14(:,i)=interp1(d14S,data14S,Depths);
    end
    if size(ind14D,1)>1
        d14D=Depths(ind14D);
        data14D=AnomD14(ind14D,i);
        AnomD14(:,i)=interp1(d14D,data14D,Depths);
    end
    if size(indD14,1)>1
        dD14=Depths(indD14);
        dataD14=D14(indD14,i);
        D14(:,i)=interp1(dD14,dataD14,Depths);
    end
    if size(ind15,1)>1
        d15=Depths(ind15);
        data15=AnomT15(ind15,i);
        AnomT15(:,i)=interp1(d15,data15,Depths);
    end
    if size(ind15S,1)>1
        d15S=Depths(ind15S);
        data15S=AnomS15(ind15S,i);
        AnomS15(:,i)=interp1(d15S,data15S,Depths);
    end
    if size(ind15D,1)>1
        d15D=Depths(ind15D);
        data15D=AnomD15(ind15D,i);
        AnomD15(:,i)=interp1(d15D,data15D,Depths);
    end
    if size(indD15,1)>1
        dD15=Depths(indD15);
        dataD15=D15(indD15,i);
        D15(:,i)=interp1(dD15,dataD15,Depths);
    end
    if size(ind16,1)>1
        d16=Depths(ind16);
        data16=AnomT16(ind16,i);
        AnomT16(:,i)=interp1(d16,data16,Depths);
    end
    if size(ind16S,1)>1
        d16S=Depths(ind16S);
        data16S=AnomS16(ind16S,i);
        AnomS16(:,i)=interp1(d16S,data16S,Depths);
    end
    if size(ind16D,1)>1
        d16D=Depths(ind16D);
        data16D=AnomD16(ind16D,i);
        AnomD16(:,i)=interp1(d16D,data16D,Depths);
    end
    if size(indD16,1)>1
        dD16=Depths(indD16);
        dataD16=D16(indD16,i);
        D16(:,i)=interp1(dD16,dataD16,Depths);
    end
    clear ind1 d1 data1 ind2 d2 data2 ind3 d3 data3 ind4 d4 data4 ind5 d5 data5...
        ind6 d6 data6 ind7 d7 data7 ind8 d8 data8 ind9 d9 data9 ind10 d10 data10...
        ind11 d11 data11 ind12 d12 data12 ind13 d13 data13 ind14 d14 data14...
        ind15 d15 data15 ind16 d16 data16 ind1S d1S data1S ind2S d2S data2S...
        ind3S d3S data3S ind4 d4S data4S ind5S d5S data5S...
        ind6S d6S data6S ind7 d7S data7S ind8S d8S data8S ind9S d9S data9S...
        ind10S d10S data10S ind11S d11S data11S ind12S d12S data12S ind13S...
        d13S data13S ind14S d14SS data14S ind15S d15S data15S ind16S d16S data16S...
        ind1D d1D data1D ind2D d2D data2D...
        ind3D d3D data3D ind4D d4D data4D ind5D d5D data5D...
        ind6D d6D data6D ind7D d7D data7D ind8D d8D data8D ind9D d9D data9D...
        ind10D d10D data10D ind11D d11D data11D ind12D d12D data12D ind13D...
        d13D data13D ind14D d14D data14D ind15D d15D data15D ind16D d16S data16D
end

%% Define figure Parameters
TMin=-3;
TMax=3;
DMin=0;
DMax=600;
DMax2=200;
SMin=-1;
SMax=1;
DenMin=-1;
DenMax=1;
DenseMin=23;
DenseMax=28;
font=24;
txtc='k';

%% Figures
%Figure 1: Temperature anomalies
fig=figure(1);
fig.Color = [1 1 1];
t=tiledlayout(4,4);
t.TileSpacing='compact';

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT1);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.XTickLabel=[];
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Winter 2013-2014';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT5);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2014-2015';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT9);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XColor='k';
ax.YColor='k';
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2015-2016';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT13);
c=colorbar;
c.Color='k';
c.FontSize=font;
c.Label.String=['TAnom (' char(176) 'C)'];
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2016-2017';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT2);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XColor='k';
ax.YColor='k';
ax.XTickLabel=[];
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Spring 2014';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT6);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2015';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT10);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2016';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT14);
c=colorbar;
c.Color='k';
c.FontSize=font;
c.Label.String=['TAnom (' char(176) 'C)'];
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2017';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT3);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Summer 2014';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT7);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2015';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT11);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2016';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT15);
c=colorbar;
c.Color='k';
c.FontSize=font;
c.Label.String=['TAnom (' char(176) 'C)'];
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2017';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT4);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 

ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2014';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT8);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2015';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT12);
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2016';
ax.Title.Color='k';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomT16);
c=colorbar;
c.Color='k';
c.FontSize=font;
c.Label.String=['TAnom (' char(176) 'C)'];
ax.Colormap=([0.8 0.8 0.8;redblue]);
ax.CLim=[TMin,TMax];
ax.YLim=[DMin,DMax]; 
ax.YTickLabel=[];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2017';
ax.Title.Color='k';
ax.TickLength=[0,0];

%Figure 2: Salinity Anomalies
fig=figure(2);
fig.Color = [1 1 1];
t=tiledlayout(4,4);
t.TileSpacing='compact';

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS1);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.XTickLabel=[];
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Winter 2013-2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS5);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2014-2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS9);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2015-2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS13);
c=colorbar;
c.FontSize=font;
c.Label.String='SAnom';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Winter 2016-2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS2);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.XTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Spring 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS6);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS10);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS14);
c=colorbar;
c.FontSize=font;
c.Label.String='SAnom';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Spring 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS3);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.XTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='Summer 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS7);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS11);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS15);
c=colorbar;
c.FontSize=font;
c.Label.String='SAnom';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='Summer 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS4);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS8);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.YTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS12);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.YTickLabel=[];
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomS16);
c=colorbar;
c.FontSize=font;
c.Label.String='Salinity Anomaly';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[SMin,SMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.YTickLabel=[];
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='Fall 2017';
ax.TickLength=[0,0];

%Figure 3: Density Anomalies
fig=figure(3);
fig.Color = [1 1 1];
t=tiledlayout(4,4);
t.TileSpacing='compact';

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD1);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='DJF 2013-2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD5);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2014-2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD9);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2015-2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD13);
c=colorbar;
c.FontSize=font;
c.Label.String='Density Anomaly';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2016-2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD2);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='MAM 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD6);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD10);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD14);
c=colorbar;
c.FontSize=font;
c.Label.String='Density Anomaly';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD3);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='JJA 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD7);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD11);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD15);
c=colorbar;
c.FontSize=font;
c.Label.String='Density Anomaly';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD4);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.XLabel.String='Longitude';
ax.Title.String='SON 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD8);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD12);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,AnomD16);
c=colorbar;
c.FontSize=font;
c.Label.String='Density Anomaly';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenMin,DenMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2017';
ax.TickLength=[0,0];

%Figure 4: Density
fig=figure(4);
fig.Color = [1 1 1];
t=tiledlayout(4,4);
t.TileSpacing='compact';

nexttile
ax=gca;
imagesc(Longitudes,Depths,D1);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='DJF 2013-2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D5);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2014-2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D9);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2015-2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D13);
c=colorbar;
c.FontSize=font;
c.Label.String='Density';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='DJF 2016-2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D2);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='MAM 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D6);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D10);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D14);
c=colorbar;
c.FontSize=font;
c.Label.String='Density';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='MAM 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D3);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.Title.String='JJA 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D7);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D11);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D15);
c=colorbar;
c.FontSize=font;
c.Label.String='Density';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.String='JJA 2017';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D4);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.XLabel.String='Longitude';
ax.Title.String='SON 2014';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D8);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2015';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D12);
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2016';
ax.TickLength=[0,0];

nexttile
ax=gca;
imagesc(Longitudes,Depths,D16);
c=colorbar;
c.FontSize=font;
c.Label.String='Density';
ax.Colormap=([0.7 0.7 0.7;customcolormap_preset('pasteljet')]);
ax.CLim=[DenseMin,DenseMax];
ax.YLim=[DMin,DMax]; 
c.Color=txtc;
ax.XColor=txtc;
ax.YColor=txtc;
ax.Title.Color=txtc;
ax.YDir='reverse';
ax.FontSize=font;
ax.LabelFontSizeMultiplier = 1.2;
ax.XLabel.String='Longitude';
ax.Title.String='SON 2017';
ax.TickLength=[0,0];


%% Contour Sections
figure(5)
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Longitudes,Depths,AnomT7,...
    [-3.5 3 -2.5 -2 -1.5 -0.5 -0.1 0 0.1 0.5 1 1.5 2 2.5 3 3.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-3 -2 -1 -0.5 -.01 0 0.1 0.5 1 2 3];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',14)
c=colorbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.YLim=[DMin,DMax];
ax.CLim=[TMin,TMax];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.Box='on';
ax.LineWidth=1.5;

figure(6)
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Longitudes,Depths,AnomT8,...
    [-3.5 3 -2.5 -2 -1.5 -0.5 -0.1 0 0.1 0.5 1 1.5 2 2.5 3 3.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-3 -2 -1 -0.5 -.01 0 0.1 0.5 1 2 3];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',14)
c=colorbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.YLim=[DMin,DMax];
ax.CLim=[TMin,TMax];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.Box='on';
ax.LineWidth=1.5;

figure(7)
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Longitudes,Depths,AnomT11,...
    [-3.5 3 -2.5 -2 -1.5 -0.5 -0.1 0 0.1 0.5 1 1.5 2 2.5 3 3.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-3 -2 -1 -0.5 -.01 0 0.1 0.5 1 2 3];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',14)
c=colorbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.YLim=[DMin,DMax];
ax.CLim=[TMin,TMax];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.Box='on';
ax.LineWidth=1.5;

figure(8)
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Longitudes,Depths,AnomT15,...
    [-3.5 3 -2.5 -2 -1.5 -0.5 -0.1 0 0.1 0.5 1 1.5 2 2.5 3 3.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-3 -2 -1 -0.5 -.01 0 0.1 0.5 1 2 3];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',14)
c=colorbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.YLim=[DMin,DMax];
ax.CLim=[TMin,TMax];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.Box='on';
ax.LineWidth=1.5;

figure(8)
fig.Color = [1 1 1];
ax=gca;
[con1,con]=contourf(Longitudes,Depths,SigmaAnomT15,...
    [-3.5 3 -2.5 -2 -1.5 -0.5 -0.1 0 0.1 0.5 1 1.5 2 2.5 3 3.5]);
con.ShowText='on';
con.LabelSpacing=220;
con.TextList=[-3 -2 -1 -0.5 -.01 0 0.1 0.5 1 2 3];
con.LineColor='k';
con.LineWidth=0.1;
clabel(con1,con,'FontSize',14)
c=colorbar;
c.Color='k';
c.FontSize=22;
c.Label.String=['Temperature Anomaly (' char(176) 'C)'];
ax.Colormap=([0 0 0;redblue]);
ax.YLim=[DMin,DMax];
ax.CLim=[TMin,TMax];
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='Depth (m)';
ax.TickLength=[0,0];
ax.Box='on';
ax.LineWidth=1.5;