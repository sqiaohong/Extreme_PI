%%
clear
clc

%%
addpath('F:\data_code\input\', ...
        'F:\data_code\output\');

load('F:\data_code\output\PI_region_GSATconstraint_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

post_all   = post_all   - mean(post_all(idx,:,:),1);
prior      = prior      - mean(prior(idx,:,:),1);

post_var  = squeeze(std(post_all,0,2));
prior_var = squeeze(std(prior,0,2));

post_best  = squeeze(post_all(:,1,:));
prior_best = squeeze(mean(prior,2));

post_95   = squeeze(prctile(post_all,95,2));
post_5    = squeeze(prctile(post_all,5,2));
prior_95  = squeeze(prctile(prior,95,2));
prior_5   = squeeze(prctile(prior,5,2));

mid = (post_5 + post_95) / 2; 
L   = (post_95 - post_5);   
L2  = 1 * L; 
post_5_new  = mid - 0.5 * L2;
post_95_new = mid + 0.5 * L2;

post_95_2024    = squeeze(post_95_new(end-1,:));
post_5_2024     = squeeze(post_5_new(end-1,:));
post_best_2024     = post_best(end-1,:);

post_uncertain90  = post_95_2024 - post_5_2024;
pvalue= post_5_2024>0;

for i=1:44
    post1=squeeze(post_all(:,2:100,i));
    pvalue=ttest(post1,0);
end


%%
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';
zones_file='F:\data_code\input\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';
[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');
S = S_all(ia);
nR = numel(S);
pvalue=pvalue(ib);

%%
load('F:\data_code\input\color20_1207.mat');
cmap = (color);
cmin = -0.1;
cmax = 0.1;

r = post_best_2024(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:); 

%%
subplot(2,2,1); hold on;

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

%%
colormap(cmap);
caxis([cmin cmax]);
title('Constrained 2024 PI forced response');
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'a','FontWeight','bold')


%%
load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);
ERA5   = series_ERA5_region(end-1,:)   - mean(series_ERA5_region(idx,:),1);
residual=ERA5-post_best_2024;

%% 
load('F:\data_code\input\color20_1207.mat');
cmap = (color);
cmin = -0.1;
cmax = 0.1;
r = residual(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:); 

%%
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

%%
subplot(2,2,2); hold on;

for i = 1:nR
    plot(S(i).Lon, S(i).Lat, 'Color',[0.6 0.6 0.6], 'LineWidth',0.7);
end

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

%% 
colormap(cmap);
caxis([cmin cmax]);
title('2024 PI anomaly relative to 1981-2010');
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'b','FontWeight','bold')



%%
clear
clc

addpath('F:\data_code\input\', ...
        'F:\data_code\output\');

load('F:\data_code\output\PI_region_GSATconstraint_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

post_all   = post_all   - mean(post_all(idx,:,:),1);
prior      = prior      - mean(prior(idx,:,:),1);

post_var  = squeeze(std(post_all,0,2));
prior_var = squeeze(std(prior,0,2));

post_best  = squeeze(post_all(:,1,:));
prior_best = squeeze(mean(prior,2));

post_95   = squeeze(prctile(post_all,95,2));
post_5    = squeeze(prctile(post_all,5,2));
prior_95  = squeeze(prctile(prior,95,2));
prior_5   = squeeze(prctile(prior,5,2));

mid = (post_5 + post_95) / 2; 
L   = (post_95 - post_5);  
L2  = 1 * L; 
post_5_new  = mid - 0.5 * L2;
post_95_new = mid + 0.5 * L2;

post_95_2024    = squeeze(post_95_new(end,:));
post_5_2024     = squeeze(post_5_new(end,:));
post_best_2024     = post_best(end,:);

post_uncertain90  = post_95_2024 - post_5_2024;
pvalue= post_5_2024>0;

for i=1:44
    post1=squeeze(post_all(:,2:100,i));
    pvalue=ttest(post1,0);
end


%%
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='F:\data_code\input\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');

S         = S_all(ia);
nR = numel(S);
pvalue=pvalue(ib);

%%
load('F:\data_code\input\color20_1207.mat');
cmap = (color);

cmin = -0.1;
cmax = 0.1;
r = post_best_2024(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:);  

%%
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

%% 
subplot(2,2,3); hold on;

for i = 1:nR
    plot(S(i).Lon, S(i).Lat, 'Color',[0.6 0.6 0.6], 'LineWidth',0.7);
end

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

%%
colormap(cmap);
caxis([cmin cmax]);
title('Constrained 2025 PI forced response');
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'c','FontWeight','bold')



%%
load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

ERA5   = series_ERA5_region(end,:)   - mean(series_ERA5_region(idx,:),1);
residual=ERA5-post_best_2024;

%% 
load('F:\data_code\input\color20_1207.mat'); 
cmap = (color);

cmin = -0.1;
cmax = 0.1;
r = residual(ib);
r(r < cmin) = cmin;
r(r > cmax) = cmax;

idx = round((r - cmin) / (cmax - cmin) * (size(cmap,1)-1)) + 1;
idx = max(1, min(idx, size(cmap,1)));
var_colors = cmap(idx,:); 

%%
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

%%
subplot(2,2,4); hold on;
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

%%
colormap(cmap);
caxis([cmin cmax]);
cb = colorbar;
title('2025 PI anomaly relative to 1981-2010');
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'d','FontWeight','bold')

