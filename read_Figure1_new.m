%%
clear
clc

%% load data
addpath('F:\data_code\input\', ...
        'F:\data_code\output\', ...
        'F:\data_code\function\');

load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

ERA5   = series_ERA5_region(end-1,:)   - mean(series_ERA5_region(idx,:),1);

%% Read Shapefile match region and data
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='F:\data_code\input\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');
S = S_all(ia);
nR = numel(S);

%% Setting Colormap
load('K:\KCC\new\CA_code_review_submitted\matlab\color20_1207.mat');
cmap = (color);
cmin = -0.1;
cmax = 0.1;

r = ERA5(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:); 

%% Plot
subplot(3,1,1); hold on;

cen_lon = nan(nR,1); 
cen_lat = nan(nR,1); 
radius  = nan(nR,1);
for i = 1:nR
    lon = S(i).Lon; 
    lat = S(i).Lat;
    ok  = ~isnan(lon) & ~isnan(lat);
    lon = lon(ok); 
    lat = lat(ok);
    try
        p = polyshape(lon,lat);
        [cx,cy] = centroid(p);
    catch
        cx = mean(lon); 
        cy = mean(lat);
    end
    cen_lon(i) = cx; 
    cen_lat(i) = cy;

    radius(i) = 0.40 * min(max(lon)-min(lon), max(lat)-min(lat));
end

for i = 1:nR
    plot(S(i).Lon, S(i).Lat, 'Color',[0.6 0.6 0.6], 'LineWidth',0.7);
end

title('Change in variability (post/prior std ratio)','FontSize',14);

theta = linspace(0, 2*pi, 80);

for i = 1:nR
    xc = cen_lon(i);
    yc = cen_lat(i);
    r0 = radius(i);

    x_circ = xc + r0 * cos(theta);
    y_circ = yc + r0 * sin(theta);

    faceColor = var_colors(i,:);

    h_circle = patch(x_circ, y_circ, 1, ...
        'FaceColor', faceColor, ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.8);

    text(xc, yc, S(i).label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',8, ...
        'FontName','Times New Roman', ...
        'FontWeight','normal', ...
        'Color','k', ...
        'BackgroundColor','none', ...
        'Margin',0.1);
end

zones_file2 = 'F:\data_code\input\regionals.shp';
region_line = shaperead(zones_file2, 'UseGeoCoords', true);
for i = 1:9
    plot(region_line(i).Lon, region_line(i).Lat, 'color',[0.3 0.3 0.3])
end

box on
set(gca,'XtickLabel',{' '},'YtickLabel',{' '})

%% Add colorbar
colormap(cmap);
caxis([cmin cmax]);
cb = colorbar('eastoutside');
title('2024 PI anomaly relative to 1981-2010');


%%
clear
clc

%% load data
addpath('F:\data_code\input\', ...
        'F:\data_code\output\', ...
        'F:\data_code\function\');

load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

ERA5   = series_ERA5_region(end,:)   - mean(series_ERA5_region(idx,:),1);

%% Read Shapefile match region and data
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='F:\data_code\input\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');
S = S_all(ia);
nR = numel(S);

%% Setting Colormap
load('K:\KCC\new\CA_code_review_submitted\matlab\color20_1207.mat');
cmap = (color);
cmin = -0.1;
cmax = 0.1;

r = ERA5(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:); 

%% Plot
subplot(3,1,2); hold on;

cen_lon = nan(nR,1); 
cen_lat = nan(nR,1); 
radius  = nan(nR,1);
for i = 1:nR
    lon = S(i).Lon; 
    lat = S(i).Lat;
    ok  = ~isnan(lon) & ~isnan(lat);
    lon = lon(ok); 
    lat = lat(ok);
    try
        p = polyshape(lon,lat);
        [cx,cy] = centroid(p);
    catch
        cx = mean(lon); 
        cy = mean(lat);
    end
    cen_lon(i) = cx; 
    cen_lat(i) = cy;

    radius(i) = 0.40 * min(max(lon)-min(lon), max(lat)-min(lat));
end

for i = 1:nR
    plot(S(i).Lon, S(i).Lat, 'Color',[0.6 0.6 0.6], 'LineWidth',0.7);
end

title('Change in variability (post/prior std ratio)','FontSize',14);

theta = linspace(0, 2*pi, 80);

for i = 1:nR
    xc = cen_lon(i);
    yc = cen_lat(i);
    r0 = radius(i);

    x_circ = xc + r0 * cos(theta);
    y_circ = yc + r0 * sin(theta);

    faceColor = var_colors(i,:);

    h_circle = patch(x_circ, y_circ, 1, ...
        'FaceColor', faceColor, ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.8);

    text(xc, yc, S(i).label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',8, ...
        'FontName','Times New Roman', ...
        'FontWeight','normal', ...
        'Color','k', ...
        'BackgroundColor','none', ...
        'Margin',0.1);
end

zones_file2 = 'F:\data_code\input\regionals.shp';
region_line = shaperead(zones_file2, 'UseGeoCoords', true);
for i = 1:9
    plot(region_line(i).Lon, region_line(i).Lat, 'color',[0.3 0.3 0.3])
end

box on
set(gca,'XtickLabel',{' '},'YtickLabel',{' '})

%% Add colorbar
colormap(cmap);
caxis([cmin cmax]);
cb = colorbar('eastoutside'); 
title('2025 PI anomaly relative to 1981-2010');



%%
clear
clc

load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
year=1950:2025;
idx=find(year>=1981 & year<=2010);
base=mean(series_ERA5_region(idx,:),1);
PI=series_ERA5_region-base;
id=PI>0;
id_num=sum(id,2);
PI1981=PI(idx,:);
id1981=PI1981>0;
id_num1981=sum(id1981,2);
id_num=sum(id,2);

subplot(3,1,3)
h1=plot(1950:2025,id_num,'k','LineWidth',1);
hold on
h2=plot([1981 2010],[mean(id_num1981) mean(id_num1981)],'--','color','r','LineWidth',1.5);
ylabel('Number of regions')
xlabel('Year');
xlim([1949 2026])
nn=roundn(mean(id_num1981),0);
text(1966,25,'1981-2010 average','Color','r')
title('(c) Number of regions with PI anomaly above zoro');

