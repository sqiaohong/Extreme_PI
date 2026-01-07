clear
clc

load('F:\data_code\output\PI_global_GSATconstraint_1950_2025.mat')

colors = [
    0.00 0.45 0.74;
    0.85 0.33 0.10;
    0.93 0.69 0.13;
    0.49 0.18 0.56;
    0.47 0.67 0.19
    ];

year_model=1950:2025;
start_ref=1981;
end_ref=2010;
idx=find(year_model>=start_ref & year_model<=end_ref);
post_all=post_all-mean(post_all(idx,:),1);
mo_da_fu=prior;
mo_da_fu=mo_da_fu-mean(mo_da_fu(idx,:),1);

p = [5 95];
q_mo_da_fu   = prctile(mo_da_fu,  p, 2);
q_post_all   = prctile(post_all,  p, 2);
uncertainty(:,:,1)=q_mo_da_fu;
uncertainty(:,:,2)=q_post_all;
un(:,1)=q_mo_da_fu(:,2)-q_mo_da_fu(:,1);
un(:,2)=q_post_all(:,2)-q_post_all(:,1);
best2(:,1)=mean(mo_da_fu,2);
best2(:,2)=post_all(:,1);

best = best2(:,2)';
p5 = uncertainty(:,1,2)';
p90 = uncertainty(:,2,2)';
t = 1950:2025;
x = [t, fliplr(t)];
y = [p5, fliplr(p90)];
fill(x, y, colors(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
    'HandleVisibility','off');
hold on
h1=plot(t, best, 'Color', colors(1,:), 'LineWidth', 1.5);

load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year=1950:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;

residual=ERA5'-best2(:,2);
for i=1:length(residual)-30
    aa=residual(i:i+30);
    std_change(i)=std(aa);
end
h2=bar(1950:2025,ERA5'-best2(:,2));
h2.FaceColor=colors(3,:);
h2.EdgeColor='none';
h2.FaceAlpha=0.4;
plot([1949 2025],[0 0], '--','Color',[0.3 0.3 0.3]);

h3=plot(1950:2025,ERA5,'.','Color',[0.5 0.5 0.5],'MarkerSize',15, 'LineWidth', 1);
h4=plot(2024,ERA5(end-1),'.','Color','r','MarkerSize',15);
h5=plot(2025,ERA5(end),'.','Color','r','MarkerSize',15);

ylim([-0.1 0.1])
xlim([1950 2026])
box on
grid off
ylabel('PI anomaly relative to 1981-2010')
legend([h3 h4 h1 h2],{'ERA5','ERA5 2024/2025','GSAT constrained PI','ERA5-GSAT constrained PI'})
set(gca,'FontName','Arial')

