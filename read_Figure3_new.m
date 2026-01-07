%%
clear
clc

%% Load data
load('F:\data_code\output\PI_global_GSATconstraint_1973_2025.mat')

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
un(:,1)=q_mo_da_fu(:,2)-q_mo_da_fu(:,1);
un(:,2)=q_post_all(:,2)-q_post_all(:,1);
best2(:,2)=mean(mo_da_fu,2);
best2(:,3)=post_all(:,1);

load('F:\data_code\input\ERA5_global_ecdf_1973_2025.mat')
year=1973:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;
best2(:,1)=ERA5;

load('F:\data_code\output\GSAT_regression_PI_region_1973_2024.mat')
yhat_ALL=yhat_ALL_global(end)

vals     = best2(end-1,:); 
err_vals = [0, un(end-1,1)/2, un(end-1,2)/2];

subplot(2,1,1)
h = bar(1,vals, 'BarWidth', 0.9); hold on

colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10;  
    0.93 0.69 0.13; 
    0.49 0.18 0.56; 
    0.47 0.67 0.19   
    ];

h(1).FaceColor = [0.3 0.3 0.3]; 
h(2).FaceColor = colors(2,:); 
h(3).FaceColor =  colors(1,:);  
h(1).FaceAlpha  =0.5;
h(2).FaceAlpha  =0.5;
h(3).FaceAlpha  =0.5;
h(1).EdgeColor='none'
h(2).EdgeColor='none'
h(3).EdgeColor='none'


xpos = zeros(1, numel(h));
for i = 1:numel(h)
    xpos(i) = h(i).XEndPoints; 
end

errorbar(xpos(2:end), vals(2:end), err_vals(2:end), ...
    'k', 'LineStyle','none', 'LineWidth',1.2);

set(gca,'XTick',1:3,'XTickLabel',{'ERA5','Prior','Post'});
ylabel('2024 PI anoamaly relative to 1981-2010')
plot(xpos(3),yhat_ALL_global(end),'.','MarkerSize',40,'Color','b')



%%
load('F:\data_code\output\PI_global_GSATconstraint_1950_2025.mat')
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
un2(:,1)=q_mo_da_fu(:,2)-q_mo_da_fu(:,1);
un2(:,2)=q_post_all(:,2)-q_post_all(:,1);
best22(:,2)=mean(mo_da_fu,2);
best22(:,3)=post_all(:,1);

load('K:\KCC\new\CA_code_review_submitted\matlab\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year=1950:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;
best22(:,1)=ERA5;

load('F:\data_code\output\GSAT_regression_PI_region_1950_2024.mat')
yhat_ALL=yhat_ALL_global(end)
vals = best22(end-1,:);
err_vals = [0, un2(end-1,1)/2, un2(end-1,2)/2];

h = bar(2,vals, 'BarWidth', 0.9); hold on
h(1).FaceColor = [0.3 0.3 0.3];  
h(2).FaceColor = colors(2,:);  
h(3).FaceColor =  colors(1,:); 
h(1).FaceAlpha  =0.5;
h(2).FaceAlpha  =0.5;
h(3).FaceAlpha  =0.5;
h(1).EdgeColor='none'
h(2).EdgeColor='none'
h(3).EdgeColor='none'

xpos = zeros(1, numel(h));
for i = 1:numel(h)
    xpos(i) = h(i).XEndPoints; 
end

errorbar(xpos(2:end), vals(2:end), err_vals(2:end), ...
    'k', 'LineStyle','none', 'LineWidth',1.2);
set(gca,'XTick',1:3,'XTickLabel',{'ERA5','Prior','Post'});
ylabel('2024 PI anoamaly relative to 1981-2010')
h4=plot(xpos(3),yhat_ALL_global(end),'.','MarkerSize',40,'Color','b')
set(gca,'Xtick',1:2,'Xticklabel',{'1973-2025' '1950-2025'})


%%
load('F:\data_code\output\PI_global_GSATconstraint_1850_2025.mat')
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
un3(:,1)=q_mo_da_fu(:,2)-q_mo_da_fu(:,1);
un3(:,2)=q_post_all(:,2)-q_post_all(:,1);
best23(:,2)=mean(mo_da_fu,2);
best23(:,3)=post_all(:,1);

load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')
year=1950:2025;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5-base;
best23(:,1)=nan.*zeros(176,1);

vals = best23(end-1,:); 
err_vals = [0, un3(end-1,1)/2, un3(end-1,2)/2];
h = bar(3,vals, 'BarWidth', 0.9); hold on
h(1).FaceColor = [0.3 0.3 0.3];  
h(2).FaceColor = colors(2,:);  
h(3).FaceColor =  colors(1,:); 
h(1).FaceAlpha  =0.5;
h(2).FaceAlpha  =0.5;
h(3).FaceAlpha  =0.5;
h(1).EdgeColor='none'
h(2).EdgeColor='none'
h(3).EdgeColor='none'

xpos = zeros(1, numel(h));
for i = 1:numel(h)
    xpos(i) = h(i).XEndPoints; 
end
errorbar(xpos(2:end), vals(2:end), err_vals(2:end), ...
    'k', 'LineStyle','none', 'LineWidth',1.2);
set(gca,'XTick',1:3,'XTickLabel',{'ERA5','Prior','Post'});
ylabel('2025 PI anoamaly relative to 1981-2010')
set(gca,'Xtick',1:3,'Xticklabel',{'1973-2025' '1950-2025' '1850-2025'})
ylim([0 0.1])
legend([h h4],{'ERA5' 'Unconstrained' 'Constrained' 'GSAT forced predicted'})
title('2024 PI forced response')
text(0.3,0.105,'a','FontWeight','bold')
set(gca,'FontName','Arial')



%%
clear
clc
colors = [
    0.00 0.45 0.74; 
    0.85 0.33 0.10;  
    0.93 0.69 0.13; 
    0.49 0.18 0.56; 
    0.47 0.67 0.19   
    ];

load('F:\data_code\output\PI_global_GSATconstraint_1973_2025_IPMtest.mat')

year = 1973:2025;
idx  = (year >= 1981 & year <= 2010);

post_all   = post_all   - mean(post_all(idx,:,:),1);
prior      = prior      - mean(prior(idx,:,:),1);
pseudo_obs = pseudo_obs - mean(pseudo_obs(idx,:),1);

post_var  = squeeze(std(post_all,0,2));
prior_var = squeeze(std(prior,0,2));

post_best  = squeeze(post_all(:,1,:));
prior_best = squeeze(mean(prior,2));

bias_post  = abs(post_best  - pseudo_obs);
bias_prior = abs(prior_best - pseudo_obs);

bias_2024(:,1) = bias_prior(end-1,:);
bias_2024(:,2) = bias_post(end-1,:);
bias_mean = squeeze(mean(bias_2024,2));

post_95=squeeze(prctile(post_all,95,2));
post_5=squeeze(prctile(post_all,5,2));
prior_95=squeeze(prctile(prior,95,2));
prior_5=squeeze(prctile(prior,5,2));


err  = squeeze(post_all(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,2,1) = sqrt(bias_2024(:,2).^2 + err_var');

err  = squeeze(prior(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,1,1) = sqrt(bias_2024(:,1).^2 + err_var');


mid = (post_5 + post_95) / 2;      
L   = (post_95 - post_5);           
L2 =   L;             
post_5_new  = mid - 0.5 * L2;
post_95_new = mid + 0.5 * L2;


pseudo_obs_2024=squeeze(pseudo_obs(end-1,:));
post_95_2024=squeeze(post_95(end-1,:));
post_5_2024=squeeze(post_5(end-1,:));

in_range = (pseudo_obs_2024 >= post_5_2024) & (pseudo_obs_2024 <= post_95_2024);
in_count =100*sum(in_range, 2)/22;
post_uncertain90=post_95_2024-post_5_2024;

prior_95_2024=squeeze(prior_95(end-1,:));
prior_5_2024=squeeze(prior_5(end-1,:));
prior_uncertain90=prior_95_2024-prior_5_2024;

uncertain90(:,1)=(prior_uncertain90);
uncertain90(:,2)=(post_uncertain90);

uncertain5(:,1)=prior_5_2024;
uncertain5(:,2)=post_5_2024;

uncertain95(:,1)=prior_95_2024;
uncertain95(:,2)=post_95_2024;

uncertain90_reduction(:,1)=100*(1-uncertain90(:,2)./uncertain90(:,1));
in_count_obs(:,1)=in_range;

%%
load('F:\data_code\output\PI_global_GSATconstraint_1950_2025_IPMtest.mat')

year = 1950:2025;
idx  = (year >= 1981 & year <= 2010);

post_all   = post_all   - mean(post_all(idx,:,:),1);
prior      = prior      - mean(prior(idx,:,:),1);
pseudo_obs = pseudo_obs - mean(pseudo_obs(idx,:),1);

post_var  = squeeze(std(post_all,0,2));
prior_var = squeeze(std(prior,0,2));

post_best  = squeeze(post_all(:,1,:));
prior_best = squeeze(mean(prior,2));

bias_post  = abs(post_best  - pseudo_obs);
bias_prior = abs(prior_best - pseudo_obs);

bias_2024(:,1) = bias_prior(end-1,:);
bias_2024(:,2) = bias_post(end-1,:);
bias_mean = squeeze(mean(bias_2024,2));

err  = squeeze(post_all(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,2,2) = sqrt(bias_2024(:,2).^2 + err_var');

err  = squeeze(prior(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,1,2) = sqrt(bias_2024(:,1).^2 + err_var');

post_95=squeeze(prctile(post_all,95,2));
post_5=squeeze(prctile(post_all,5,2));
prior_95=squeeze(prctile(prior,95,2));
prior_5=squeeze(prctile(prior,5,2));


mid = (post_5 + post_95) / 2;  
L   = (post_95 - post_5); 
L2 =   L;         
post_5_new  = mid - 0.5 * L2;
post_95_new = mid + 0.5 * L2;


pseudo_obs_2024=squeeze(pseudo_obs(end-1,:));
post_95_2024=squeeze(post_95(end-1,:));
post_5_2024=squeeze(post_5(end-1,:));

in_range = (pseudo_obs_2024 >= post_5_2024) & (pseudo_obs_2024 <= post_95_2024);
in_count =100*sum(in_range, 2)/22;
post_uncertain90=post_95_2024-post_5_2024;

prior_95_2024=squeeze(prior_95(end-1,:));
prior_5_2024=squeeze(prior_5(end-1,:));
prior_uncertain90=prior_95_2024-prior_5_2024;

uncertain90(:,1)=(prior_uncertain90);
uncertain90(:,2)=(post_uncertain90);

uncertain5(:,1)=prior_5_2024;
uncertain5(:,2)=post_5_2024;

uncertain95(:,1)=prior_95_2024;
uncertain95(:,2)=post_95_2024;
uncertain90_reduction(:,2)=100*(1-uncertain90(:,2)./uncertain90(:,1));
in_count_obs(:,2)=in_range;


%%
load('F:\data_code\output\PI_global_GSATconstraint_1850_2025_IPMtest.mat')
year = 1850:2025;
idx  = (year >= 1981 & year <= 2010);

post_all   = post_all   - mean(post_all(idx,:,:),1);
prior      = prior      - mean(prior(idx,:,:),1);
pseudo_obs = pseudo_obs - mean(pseudo_obs(idx,:),1);

post_var  = squeeze(std(post_all,0,2));
prior_var = squeeze(std(prior,0,2));

post_best  = squeeze(post_all(:,1,:));
prior_best = squeeze(mean(prior,2));

bias_post  = abs(post_best  - pseudo_obs);
bias_prior = abs(prior_best - pseudo_obs);

bias_2024(:,1) = bias_prior(end-1,:);
bias_2024(:,2) = bias_post(end-1,:);
bias_mean = squeeze(mean(bias_2024,2));

err  = squeeze(post_all(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,2,3) = sqrt(bias_2024(:,2).^2 + err_var');

err  = squeeze(prior(end-1,:,:))- pseudo_obs(end-1,:);      % error of each model
err_var = var(err, 0,1); 
rmse(:,1,3) = sqrt(bias_2024(:,1).^2 + err_var');

post_95=squeeze(prctile(post_all,95,2));
post_5=squeeze(prctile(post_all,5,2));
prior_95=squeeze(prctile(prior,95,2));
prior_5=squeeze(prctile(prior,5,2));

mid = (post_5 + post_95) / 2;  
L   = (post_95 - post_5); 
L2 =   1.3*L;         
post_5_new  = mid - 0.5 * L2;
post_95_new = mid + 0.5 * L2;

pseudo_obs_2024=squeeze(pseudo_obs(end-1,:));
post_95_2024=squeeze(post_95_new(end-1,:));
post_5_2024=squeeze(post_5_new(end-1,:));

in_range = (pseudo_obs_2024 >= post_5_2024) & (pseudo_obs_2024 <= post_95_2024);
in_count =100*sum(in_range, 2)/22;
post_uncertain90=post_95_2024-post_5_2024;

prior_95_2024=squeeze(prior_95(end-1,:));
prior_5_2024=squeeze(prior_5(end-1,:));
prior_uncertain90=prior_95_2024-prior_5_2024;

uncertain90(:,1)=(prior_uncertain90);
uncertain90(:,2)=(post_uncertain90);

uncertain5(:,1)=prior_5_2024;
uncertain5(:,2)=post_5_2024;

uncertain95(:,1)=prior_95_2024;
uncertain95(:,2)=post_95_2024;

uncertain90_reduction(:,3)=100*(1-uncertain90(:,2)./uncertain90(:,1));
in_count_obs(:,3)=in_range;

subplot(2,1,2)
h=bar(uncertain90_reduction,'BarWidth', 0.9)
set(gca,'Xtick',1:22, 'Xticklabel', common_models);
h(1).FaceColor=colors(3,:);
h(2).FaceColor=colors(4,:);
h(3).FaceColor=colors(5,:);
h(1).FaceAlpha  =0.5;
h(2).FaceAlpha  =0.5;
h(3).FaceAlpha  =0.5;
h(1).EdgeColor='none'
h(2).EdgeColor='none'
h(3).EdgeColor='none'

ylabel('Reduction in uncertainty range (%)')

sig_idx = in_count_obs;  
sig_idx = sig_idx ~= 0;              

data_all = uncertain90_reduction;    
data_hatch = nan(size(data_all));     
data_hatch(sig_idx) = data_all(sig_idx);

hold on
h_sig = bar(data_hatch,'BarWidth', 0.9); 

set(h_sig(1),'FaceColor','none','EdgeColor','none');
set(h_sig(2),'FaceColor','none','EdgeColor','none');
set(h_sig(3),'FaceColor','none','EdgeColor','none');

if any(sig_idx(:,1))
    bb=hatchfill2(h_sig(1),'single', ...
        'HatchAngle',45, ...
        'HatchDensity',50, ...
        'HatchColor',[0.3 0.3 0.3]);
end

if any(sig_idx(:,2))
    hatchfill2(h_sig(2),'single', ...
        'HatchAngle',45, ...
        'HatchDensity',50, ...
        'HatchColor',[0.3 0.3 0.3]);
end

if any(sig_idx(:,3))
    hatchfill2(h_sig(3),'single', ...
        'HatchAngle',45, ...
        'HatchDensity',50, ...
        'HatchColor',[0.3 0.3 0.3]);
end

RMSE_mean=squeeze(mean(rmse,1));
RMSE_ratio=RMSE_mean(2,:)./RMSE_mean(1,:);

legend([h(1) h(2) h(3) bb],{'1973-2025' '1950-2025'  'Adjust 1850-2025' 'Pseudo in constrained 90% confident coverage '})
title('Imperfect model test for KCC constraint')
ylim([0 40])

RMSE_ratio=roundn(RMSE_ratio,-2);
text(11,38,num2str(RMSE_ratio(1)),'color',colors(3,:))
text(11,36,num2str(RMSE_ratio(2)),'color',colors(4,:))
text(11,34,num2str(RMSE_ratio(3)),'color',colors(5,:))

text(-0.2,41,'b','FontWeight','bold')
set(gca,'FontName','Arial')




