%% 0.25dg
clear
clc

addpath('F:\data_code\input\', 'F:\data_code\output\', 'F:\data_code\function\');

load('F:\data_code\input\ERA5_Rx1day_ecdf_1950_2025.mat')
cdf_vals=ERA5_ecdf;
cdf_vals(cdf_vals<0)=nan; 
load('F:\data_code\input\global_mask_ERA5.mat')
load('F:\data_code\input\ERA5_lon_lat_1950.mat')
mask1=flipud(mask');
mask2(:,1:720)=mask1(:,721:end);
mask2(:,721:1440)=mask1(:,1:720);
mask2(mask2==0)=nan;
id=find(lat<-60);
mask3=mask2;
mask3(id,:)=nan;
load('F:\data_code\input\ERA5_area_weight.mat');
for i=1:75
    data1=cdf_vals(:,:,i).*mask3';
    area_weight3=area_weight.*mask3';
    data1=data1.*area_weight3;
    series_ERA5(i)=nansum(data1(:))./nansum(area_weight3(:));
end

year=1950:2024;
id1=find(year==1981);
id2=find(year==2010);
base=mean(series_ERA5(id1:id2));
cdf2024=series_ERA5(end)-base;

ERA5=load('F:\data_code\input\ERA5_Rx1day_1950_2025.mat');

year=1950:2024;
id1=find(year==1981);
id2=find(year==2010);

ERA5_cdf=NaN.*zeros(size(ERA5.Rx1day,1),size(ERA5.Rx1day,2),size(ERA5.Rx1day,3));
cdf_base=nanmean(cdf_vals(:,:,id1:id2),3);

cdf2024=0.094;
for i=1:size(ERA5.Rx1day,1)
    i
    for j=1:size(ERA5.Rx1day,2)
        rx=squeeze(ERA5.Rx1day(i,j,id1:id2));
        rx1=squeeze(ERA5.Rx1day(i,j,:));
        if isnan(rx)~=1
            base1=cdf_base(i,j);
            cdf20241=cdf2024+base1;
            [ERA5_value1981(i,j,:), P_sorted] = ecdf_plot_inverse(rx1, base1);
            [ERA5_value2024(i,j,:), P_sorted] = ecdf_plot_inverse(rx1,cdf20241);
        end
    end
end

load('E:\PCIC\grid\ERA5_area_weight.mat');
change=100*(ERA5_value2024-ERA5_value1981)./ERA5_value1981;
% change1=change.*mask3'.*area_weight3;
% global_change=nansum(change1(:)) ./ nansum(area_weight3(:)); 
median_value=nanmedian(change(:));





%% 2.5dg
clear
clc

addpath('F:\data_code\input\', 'F:\data_code\output\', 'F:\data_code\function\');

load('F:\data_code\input\ERA5_Rx1day_ecdf_1950_2025_25dg.mat')
lat_bins = -90:2.5:90;
lon_bins = -180:2.5:177.5;
ERA5_cdf(ERA5_cdf==0)=nan;
load('F:\data_code\input\ERA5_area_weight_25dg.mat');
load('F:\data_code\input\global_mask_ERA5_25dg.mat')
mask_global=mask;
mask_global(mask_global==0)=nan;
id=find(lat_bins<-60);
mask_global(:,id)=nan;

area_weight=area_weight.*mask_global;
for i=1:size(ERA5_cdf,3)
    data=ERA5_cdf(:,:,i)';
    data=data'.*mask_global.*area_weight;
    idxx=find(isnan(data)==0);
    series_ERA5(i)=nansum(data(:)) ./ nansum(area_weight(:));
end

year=1950:2024;
id1=find(year==1981);
id2=find(year==2010);
base=mean(series_ERA5(id1:id2));
cdf2024=series_ERA5(end-1)-base;

ERA5=load('F:\data_code\input\ERA5_Rx1day_1950_2025.mat');


load('F:\data_code\input\ERA5_Rx1day_1950_2025.mat')
lat_bins = -90:2.5:90;
lon_bins = -180:2.5:177.5;

[glon, glat] = meshgrid(lon_bins, lat_bins);
glon = glon';
glat = glat';

load('F:\data_code\input\ERA5_lon_lat_1950.mat')
[glon1, glat1] = meshgrat(lon, lat);

data=Rx1day;
data= reshape(data,1440*721,76);
ERA5.lon=reshape(glon1,1440*721,1);
ERA5.lat=reshape(glat1,1440*721,1);


ERA5_Rx1day = zeros(144, 73, 76) * nan;

glon_edges = [glon(1:end-1, 1); glon(end, 1)]; 
glat_edges = [glat(1, 1:end-1), glat(1, end)]; 

[~, ~, lon_bin] = histcounts(ERA5.lon, glon_edges);
[~, ~, lat_bin] = histcounts(ERA5.lat, glat_edges);

valid_idx = (lon_bin > 0) & (lat_bin > 0) & (lon_bin <= 144) & (lat_bin <= 73);

for t = 1:76
    ERA5_Rx1day(:, :, t) = accumarray([lon_bin(valid_idx), lat_bin(valid_idx)], ...
        data(valid_idx, t), [144, 73], @nanmean);
end

year=1950:2024;
id1=find(year==1981);
id2=find(year==2010);

cdf_base=nanmean(ERA5_cdf(:,:,id1:id2),3);

for i=1:size(cdf_base,1)
    i
    for j=1:size(cdf_base,2)
        rx=squeeze(ERA5_Rx1day(i,j,id1:id2));
        rx1=squeeze(ERA5_Rx1day(i,j,:));
        if isnan(rx)~=1
            base1=cdf_base(i,j);
            cdf20241=cdf2024+base1;
            [ERA5_value1981(i,j,:), P_sorted] = ecdf_plot_inverse(rx1, base1);
            [ERA5_value2024(i,j,:), P_sorted] = ecdf_plot_inverse(rx1,cdf20241);
        end
    end
end


change=100*(ERA5_value2024-ERA5_value1981)./ERA5_value1981;
change1=change.*mask_global.*area_weight;
global_change=nansum(change1(:)) ./ nansum(area_weight(:)); 
median_value=nanmedian(change(:));





