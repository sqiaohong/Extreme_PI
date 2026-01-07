clear
clc

%% 1973-2024
load('F:\data_code\input\ERA5_Rx1day_ecdf_1973_2025.mat')

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
for i=1:53
    data1=cdf_vals(:,:,i).*mask3';
    idx=find(data1<=0.9);
    idx2=find(data1>0.9);
    data2=data1;
    data2(idx)=nan;
    area_weight3=area_weight.*mask3';
    area9(i)=nansum(area_weight3(idx2))./nansum(area_weight3(:));
end
area9=area9.*100;

for i=1:53
    data1=cdf_vals(:,:,i).*mask3';
    area_weight3=area_weight.*mask3';
    data1=data1.*area_weight3;
    series_ERA5(i)=nansum(data1(:))./nansum(area_weight3(:));
end

colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10; 
    0.93 0.69 0.13; 
    0.49 0.18 0.56; 
    0.47 0.67 0.19 
    ];

year=1973:2025;
subplot(2,1,1)
yyaxis left;
h1=plot(year,series_ERA5,'color',colors(1,:),'LineWidth',1,'LineStyle','-','Marker','o');
hold on
ylim([0.4 0.6])
xlim([1949 2026])
ylabel('Rx1day Probability index')
text(1947,0.61,'a','FontWeight','bold')


yyaxis right;
harea9=plot(year,area9,'-','color',colors(2,:),'LineWidth',1);
harea91=plot([1949 2026],[10 10],'--','color',[0.7 0.7 0.7],'LineWidth',1);
ylim([0 20])
legend([h1 harea9],{'PI' 'Area with PI>0.9'});
ylabel('Percentage of area (%)')
title('1973-2025')




%% 1950-2024
clear
clc

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
for i=1:76
    data1=cdf_vals(:,:,i).*mask3';
    idx=find(data1<=0.9);
    idx2=find(data1>0.9);
    data2=data1;
    data2(idx)=nan;
    area_weight3=area_weight.*mask3';
    area9(i)=nansum(area_weight3(idx2))./nansum(area_weight3(:));
end
area9=area9.*100;

for i=1:76
    data1=cdf_vals(:,:,i).*mask3';
    area_weight3=area_weight.*mask3';
    data1=data1.*area_weight3;
    series_ERA5(i)=nansum(data1(:))./nansum(area_weight3(:));
end

colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10; 
    0.93 0.69 0.13; 
    0.49 0.18 0.56; 
    0.47 0.67 0.19 
    ];
year=1950:2025;
subplot(2,1,2)
yyaxis left;
h1=plot(year,series_ERA5,'color',colors(1,:),'LineWidth',1,'LineStyle','-','Marker','o');
hold on
ylim([0.4 0.6])
xlim([1949 2026])
ylabel('Rx1day Probability index')
text(1947,0.61,'b','FontWeight','bold')

yyaxis right;
harea9=plot(year,area9,'-','color',colors(2,:),'LineWidth',1);
harea91=plot([1949 2026],[10 10],'--','color',[0.7 0.7 0.7],'LineWidth',1);
ylim([0 20])
ylabel('Percentage of area (%)')
title('1950-2025')

set(gca,'FontName','Arial')


