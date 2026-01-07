%% Setting Colormap (one unique color per IPCC region; NOT value-based)
load('F:\data_code\input\color20_1207.mat'); % provides "color"
cmap0 = color;

% build a 44-color map by interpolating your base colormap
nR = 44;
xi = linspace(1, size(cmap0,1), nR);
cmap_regions = interp1(1:size(cmap0,1), cmap0, xi, 'linear');

%% Plot: IPCC regions only (no value-based colors)
hold on;

% take only 44 IPCC land regions
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='F:\data_code\input\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');
S = S_all(ia);
nR = numel(S);

% compute centroids + a radius for label background circle (optional)
cen_lon = nan(nR,1);
cen_lat = nan(nR,1);
radius  = nan(nR,1);

for i = 1:nR
    lon = S(i).Lon;
    lat = S(i).Lat;

    ok = ~isnan(lon) & ~isnan(lat);
    lon = lon(ok);
    lat = lat(ok);

    % some shapes may include separators; polyshape can fail
    try
        p = polyshape(lon, lat);
        [cx, cy] = centroid(p);
    catch
        cx = mean(lon);
        cy = mean(lat);
    end
    cen_lon(i) = cx;
    cen_lat(i) = cy;

    radius(i) = 0.40 * min(max(lon)-min(lon), max(lat)-min(lat));
end

% draw filled IPCC polygons (each region has one fixed color)
for i = 1:nR
    faceColor = 'r';

    patch(S(i).Lon, S(i).Lat, 1, ...
        'FaceColor', [0.5 0.5 1], ...
        'EdgeColor', 'k', ...
        'LineWidth', 0.8);
end

% add region acronym labels (use Acronym field)
for i = 1:nR
    text(cen_lon(i), cen_lat(i), S(i).label, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',10, ...
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

title('IPCC reference regions','FontSize',14);

box on
set(gca,'XtickLabel',{' '},'YtickLabel',{' '})

% IMPORTANT: no value-based colormap/caxis/colorbar
% (so remove / comment out your colormap(cmap); caxis(...); colorbar(...)
