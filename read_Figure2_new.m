clear
clc

load('F:\data_code\output\PI_global_GSATconstraint_1973_2025.mat')

colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10; 
    0.93 0.69 0.13; 
    0.49 0.18 0.56; 
    0.47 0.67 0.19 
    ];

year_model=1973:2025;
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
subplot(3,1,1)

for ff = 1:2

    best = best2(:,ff)';
    p5 = uncertainty(:,1,ff)';
    p90 = uncertainty(:,2,ff)';
    t = 1973:2025;
    x = [t, fliplr(t)];
    y = [p5, fliplr(p90)];
    fill(x, y, colors(3-ff,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'HandleVisibility','off');
    hold on
    h(ff)=plot(t, best, 'Color', colors(3-ff,:), 'LineWidth', 1.5);

end

load('F:\data_code\input\ERA5_global_ecdf_1973_2025.mat')
year=1973:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;
plot([1949 2025],[0 0], '--','Color',[0.3 0.3 0.3]);
h1=plot(1973:2025,ERA5,'.','Color',[0.5 0.5 0.5],'MarkerSize',15, 'LineWidth', 1);
h2=plot(2024,ERA5(end-1),'.','Color','r','MarkerSize',15);
h3=plot(2025,ERA5(end),'.','Color','r','MarkerSize',15);

load('F:\data_code\output\GSAT_regression_PI_region_1973_2024.mat')
h4=plot(1973:2024,yhat_ALL_global(:),'-','Color',colors(5,:), 'LineWidth', 1.5);
ylim([-0.1 0.1])
xlim([1973 2025])
title('1973-2025')
box on
grid off
ylabel('PI anomaly relative to 1981-2010')
legend([h1 h2 h(1) h(2) h4],{'ERA5','ERA5 2024/2025','Unconstrained PI','GSAT constrained PI','GSAT forced predicted'})
text(1970,0.12,'a','FontWeight','bold')
set(gca,'FontName','Arial')


%%
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
subplot(3,1,2)

for ff = 1:2

    best = best2(:,ff)';
    p5 = uncertainty(:,1,ff)';
    p90 = uncertainty(:,2,ff)';
    t = 1950:2025;
    x = [t, fliplr(t)];
    y = [p5, fliplr(p90)];
    fill(x, y, colors(3-ff,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'HandleVisibility','off');
    hold on
    h(ff)=plot(t, best, 'Color', colors(3-ff,:), 'LineWidth', 1.5);

end

load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year=1950:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;

plot([1949 2025],[0 0], '--','Color',[0.3 0.3 0.3]);
h1=plot(1950:2025,ERA5,'.','Color',[0.5 0.5 0.5],'MarkerSize',15, 'LineWidth', 1);
h2=plot(2024,ERA5(end-1),'.','Color','r','MarkerSize',15);
h2=plot(2025,ERA5(end),'.','Color','r','MarkerSize',15);

load('F:\data_code\output\GSAT_regression_PI_region_1950_2024.mat')
h4=plot(1950:2024,yhat_ALL_global(:),'-','Color',colors(5,:), 'LineWidth', 1.5);
ylim([-0.1 0.1])
xlim([1949 2025])
title('1950-2025')

box on
grid off
ylabel('PI anomaly relative to 1981-2010')

text(1947,0.12,'b','FontWeight','bold')
set(gca,'FontName','Arial')







%%
clear
clc

load('F:\data_code\output\PI_global_GSATconstraint_1850_2025.mat')

colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10;  
    0.93 0.69 0.13;   
    0.49 0.18 0.56;  
    0.47 0.67 0.19   
    ];

year_model=1850:2025;
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
subplot(3,1,3)

for ff = 1:2

    best = best2(:,ff)';
    p5 = uncertainty(:,1,ff)';
    p90 = uncertainty(:,2,ff)';
    t = 1850:2025;
    x = [t, fliplr(t)];
    y = [p5, fliplr(p90)];
    fill(x, y, colors(3-ff,:), 'EdgeColor', 'none', 'FaceAlpha', 0.2, ...
        'HandleVisibility','off');
    hold on
    h(ff)=plot(t, best, 'Color', colors(3-ff,:), 'LineWidth', 1.5);

end

load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year=1950:2025;
idd1=find(year==1981);
idd2=find(year==2010);

base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;

plot([1849 2025],[0 0], '--','Color',[0.3 0.3 0.3]);

h1=plot(1950:2025,ERA5,'.','Color',[0.5 0.5 0.5],'MarkerSize',15, 'LineWidth', 1);
h2=plot(2024,ERA5(end-1),'.','Color','r','MarkerSize',15);
h2=plot(2025,ERA5(end),'.','Color','r','MarkerSize',15);
load('F:\data_code\output\GSAT_regression_PI_region_1950_2024.mat')
ylim([-0.1 0.1])
xlim([1849 2025])
title('1850-2025')
box on
grid off
ylabel('PI anomaly relative to 1981-2010')
text(1847,0.12,'c','FontWeight','bold')
set(gca,'FontName','Arial')


