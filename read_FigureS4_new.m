%%
clear
clc

load('F:\data_code\input\NAO_month_1850_2025.mat')
year_NAO=1850:2025;
idx=find(year_NAO>=1951 & year_NAO<=2024);
NAO=mean(NAO(:,idx),1);

load('F:\data_code\input\nino34_month_1950_2025.mat')
year_nino34=1950:2025;
nino34(nino34<-100)=nan;
idx=find(year_nino34>=1950 & year_nino34<=2024);
nino341=cat(1,nino34(12,1:end-1),nino34(1:2,2:end));
nino34=mean(nino341,1);


load('F:\data_code\input\PDO_month_1854_2025.mat')
year_PDO=1854:2025;
idx=find(year_PDO>=1951 & year_PDO<=2024);
PDO=mean(PDO(:,idx),1);

load('F:\data_code\input\AO_month_1950_2025.mat')
year_AO=1950:2025;
idx=find(year_AO>=1951 & year_AO<=2024);
AO=mean(AO(:,idx),1);

load('F:\data_code\input\ONI_month_1950_2025.mat')
year_ONI=1950:2025;
idx=find(year_ONI>=1951 & year_ONI<=2024);
% ONI=mean(ONI(:,idx),1);
ONI1=cat(1,ONI(12,1:end-1),ONI(1:2,2:end));
ONI=mean(ONI1,1);


load('F:\data_code\input\IOD_month_1854_2024.mat')
year_IOD=1854:2024;
idx=find(year_IOD>=1951 & year_IOD<=2024);
IOD=mean(IOD(:,idx),1);

load('F:\data_code\input\AMO_month_1854_2025.mat')
year_AMO=1854:2025;
idx=find(year_AMO>=1951 & year_AMO<=2024);
for i=1:12
    AMO2(i,:)=detrend(AMO(i,idx));
end
AMO=mean(AMO2(:,:),1);


load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year_PI=1950:2025;
idx=find(year_PI>=1951 & year_PI<=2024);
PI=series_ERA5(idx);
PI_region=series_ERA5_region(idx,:);

load('F:\data_code\output\PI_region_GSATconstraint_1950_2025.mat')
PI_forced=squeeze(post_all(idx,1,:));

year=1951:2024;
idx=find(year_PI>=1981 & year_PI<=2010);
PI=PI-mean(PI(idx));
PI_region=PI_region-mean(PI_region(idx,:),1);
PI_forced=PI_forced-mean(PI_forced(idx,:),1);

index(:,1)=(nino34(1:end-1));
index(:,2)=(PDO(1:end))-mean(PDO(idx));
index(:,3)=(IOD(1:end))-mean(IOD(idx));
index(:,4)=(NAO(1:end))-mean(NAO(idx));
index(:,5)=(AO(1:end))-mean(AO(idx));
index(:,6)=(AMO(1:end))-mean(AMO(idx));

year=1951:2024;
id=find(year>=1951 & year<=2024);
for m=1:44
    for i=1:6
        PII=PI_region(:,m)-PI_forced(:,m);
%         [corr_index_region(m,i) pvalue_region(m,i)]=corr(index(:,i),PI_region(:,m));
        [corr_index_region2(m,i) pvalue_region2(m,i)]=corr((index(:,i)),PII);
    end
end


ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='K:\KCC\new\CA_code_review_submitted\matlab\zones_layers_shp\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');

S = S_all(ia);
bias_ratio=corr_index_region2(ib,1);
rmse_ratio=corr_index_region2(ib,3);
var_ratio =corr_index_region2(ib,2);

p_bias_sel=pvalue_region2(ib,1);
p_rmse_sel=pvalue_region2(ib,3);
p_var_sel =pvalue_region2(ib,2);

nR = numel(S);


ratio_mat = [bias_ratio rmse_ratio var_ratio];
p_mat     = [p_bias_sel p_rmse_sel p_var_sel];

load('K:\KCC\new\CA_code_review_submitted\matlab\color20_1207.mat');
cmap = (color);;

cmin=-1; cmax=1;

colors_map = nan(nR,3,3);

for k=1:3
    r = ratio_mat(:,k);
    r(r<cmin)=cmin;
    r(r>cmax)=cmax;

    idx=round((r-cmin)/(cmax-cmin)*(size(cmap,1)-1))+1;
    idx=max(1,min(idx,size(cmap,1)));

    colors_map(:,k,:) = cmap(idx,:);
end


cen_lon=nan(nR,1); cen_lat=nan(nR,1); radius=nan(nR,1);

for i=1:nR
    lon=S(i).Lon; lat=S(i).Lat;
    ok=~isnan(lon)&~isnan(lat);
    lon=lon(ok); lat=lat(ok);

    try
        p=polyshape(lon,lat);
        [cx,cy]=centroid(p);
    catch
        cx=mean(lon); cy=mean(lat);
    end

    cen_lon(i)=cx; cen_lat(i)=cy;

    radius(i)=0.40 * min(max(lon)-min(lon), max(lat)-min(lat)); 
end


subplot(2,1,1); hold on;

for i=1:nR
    plot(S(i).Lon, S(i).Lat,'Color',[0.6 0.6 0.6],'LineWidth',0.7);
end

ang_edges=[0,2*pi/3,4*pi/3,2*pi];

for i=1:nR
    xc=cen_lon(i); yc=cen_lat(i); r0=radius(i);

    for k=1:3
        faceColor=squeeze(colors_map(i,k,:))';

        ang = linspace(ang_edges(k),ang_edges(k+1),40);
        x_c = xc + [0 r0*cos(ang)];
        y_c = yc + [0 r0*sin(ang)];

        h_wedge = patch(x_c, y_c, 1, ...
            'FaceColor',faceColor, ...
            'EdgeColor',[1 1 1], ...
            'LineWidth',0.8);

        if p_mat(i,k)<0.1
            hatchfill2(h_wedge,'single', ...
                'HatchAngle',45, ...
                'HatchDensity',220, ...
                'HatchColor',[0 0 0]);
        end
    end

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

zones_file = 'K:\KCC\new\CA_code_review_submitted\matlab\zones_layers_shp\regionals.shp';
region_line=shaperead(zones_file, 'UseGeoCoords', true);
for i=1:9
    plot(region_line(i).Lon,region_line(i).Lat,'color',[0.3 0.3 0.3])
end

colormap(cmap);
caxis([cmin cmax]);

cb.FontSize=12;
title('Correlation with Nino34, PDO, and IOD');
box on
set(gca,'XtickLabel',{' '},'YtickLabel',{' '})
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'a','FontWeight','bold')



%%
ipcc_file='F:\data_code\input\IPCC-WGI-reference-regions-v4.shp';
ipcc = shaperead(ipcc_file,'UseGeoCoords',true);
ipcc_name={ipcc(1:44).Acronym}';

zones_file='K:\KCC\new\CA_code_review_submitted\matlab\zones_layers_shp\zones.shp';
S_all = shaperead(zones_file,'UseGeoCoords',true);
regionNames = {S_all.label}';

[A, ia, ib] = intersect(regionNames, ipcc_name,'stable');

S = S_all(ia);
bias_ratio=corr_index_region2(ib,4);
rmse_ratio=corr_index_region2(ib,6);
var_ratio =corr_index_region2(ib,5);

p_bias_sel=pvalue_region2(ib,4);
p_rmse_sel=pvalue_region2(ib,6);
p_var_sel =pvalue_region2(ib,5);

nR = numel(S);


ratio_mat = [bias_ratio rmse_ratio var_ratio];
p_mat     = [p_bias_sel p_rmse_sel p_var_sel];

load('K:\KCC\new\CA_code_review_submitted\matlab\color20_1207.mat'); 
cmap = (color);;

cmin=-1; cmax=1;

colors_map = nan(nR,3,3);

for k=1:3
    r = ratio_mat(:,k);
    r(r<cmin)=cmin;
    r(r>cmax)=cmax;

    idx=round((r-cmin)/(cmax-cmin)*(size(cmap,1)-1))+1;
    idx=max(1,min(idx,size(cmap,1)));

    colors_map(:,k,:) = cmap(idx,:);
end


cen_lon=nan(nR,1); cen_lat=nan(nR,1); radius=nan(nR,1);

for i=1:nR
    lon=S(i).Lon; lat=S(i).Lat;
    ok=~isnan(lon)&~isnan(lat);
    lon=lon(ok); lat=lat(ok);

    try
        p=polyshape(lon,lat);
        [cx,cy]=centroid(p);
    catch
        cx=mean(lon); cy=mean(lat);
    end

    cen_lon(i)=cx; cen_lat(i)=cy;

    radius(i)=0.40 * min(max(lon)-min(lon), max(lat)-min(lat));  
end


subplot(2,1,2); hold on;

for i=1:nR
    plot(S(i).Lon, S(i).Lat,'Color',[0.6 0.6 0.6],'LineWidth',0.7);
end

ang_edges=[0,2*pi/3,4*pi/3,2*pi];

for i=1:nR
    xc=cen_lon(i); yc=cen_lat(i); r0=radius(i);

    for k=1:3
        faceColor=squeeze(colors_map(i,k,:))';

        ang = linspace(ang_edges(k),ang_edges(k+1),40);
        x_c = xc + [0 r0*cos(ang)];
        y_c = yc + [0 r0*sin(ang)];


        h_wedge = patch(x_c, y_c, 1, ...
            'FaceColor',faceColor, ...
            'EdgeColor',[1 1 1], ...
            'LineWidth',0.8);


        if p_mat(i,k)<0.1
            hatchfill2(h_wedge,'single', ...
                'HatchAngle',45, ...
                'HatchDensity',220, ...
                'HatchColor',[0 0 0]);
        end
    end


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

zones_file = 'K:\KCC\new\CA_code_review_submitted\matlab\zones_layers_shp\regionals.shp';
region_line=shaperead(zones_file, 'UseGeoCoords', true);
for i=1:9
    plot(region_line(i).Lon,region_line(i).Lat,'color',[0.3 0.3 0.3])
end

colormap(cmap);
caxis([cmin cmax]);

title('Correlation with NAO, AO, and AMO');
box on
set(gca,'XtickLabel',{' '},'YtickLabel',{' '})
set(gca,'FontName','Arial')
text(-3000000-100,1100000,'b','FontWeight','bold')



