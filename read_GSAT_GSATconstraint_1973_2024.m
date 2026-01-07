%%
clear
clc

addpath('F:\data_code\input\', 'F:\data_code\output\', 'F:\data_code\function\');
start_ref=1973;
end_ref=2024;

start_year=1973;
end_year=2024;

%% GSAT Raw models from 1850-2100
load('F:\data_code\input\tas_global_GSAT_allSSP_2p5deg.mat')
model_fut=GAST_out.model{2};
model_names = cellfun(@(x) regexp(x, 'tas_annual_average_(.+?)_', 'tokens', 'once'), ...
    model_fut, 'UniformOutput', false);
model_GSAT_name = cellfun(@(x) x{1}, model_names, 'UniformOutput', false);
GSAT_model=squeeze(GAST_out.GAST_model(:,:,2)); %% ssp245
year_model=GAST_out.year;
idx=find(year_model>=start_year & year_model<=end_year);
year_model=start_year:end_year;
for i=1:size(GSAT_model,2)
    data=GSAT_model(idx,i);
    smooth1 = ncs_random_smooth(year_model, data, 'NRepeat', 10, 'SpanYears', [5 15], 'ReturnAll', true);
    GSAT_model_smooth(:,i)=smooth1.y_smooth;
end
iy=find(year_model>=start_ref & year_model<=end_ref);
mod_ref = GSAT_model_smooth - nanmean(GSAT_model_smooth(iy,:), 1);


%% Large ensemble here for estimate obs IV Cov
load('F:\data_code\input\tas_global_GSAT_large_2p5.mat')
year_obs=1850:2024;
idx=find(year_obs>=start_year & year_obs<=end_year);
year_obs=year_obs(idx);
ln_GSAT_data=GLOBAL_large.tas_all(idx,:);
model_name_large=GLOBAL_large.model_name;
iy=find(year_obs>=start_ref & year_obs<=end_ref);
ln_GSAT_data_ref  = ln_GSAT_data  - (mean(ln_GSAT_data(iy,:),  1, 'omitnan'));


%% observation, 200 runs HadCRUT5
load('F:\data_code\input\HadCRUT5_GSAT_2p5_global.mat')
obs_200runs =HadCRUT5_GSAT_2p5.tas_global(1:12*length(1850:2024),:);
obs_200runs = reshape(obs_200runs, 12, length(1850:2024), 200);
obs_200runs=squeeze(mean(obs_200runs,1));
obs_200runs=obs_200runs(idx,:);
obs_200runs_ref= obs_200runs- (mean(obs_200runs(iy,:),1, 'omitnan'));
obs_ref =squeeze(mean(obs_200runs_ref,2));

%%
scale_m=1;
obs_best          = obs_ref(:);        
runs_data(:,1,:)  = obs_200runs_ref;  
ln_data(:,1,:)    = ln_GSAT_data_ref;
mod_his(:,1,:)    = mod_ref;     
runs_data = reshape(runs_data, size(runs_data,1)*size(runs_data,2), size(runs_data,3));
obs_mea_cov = cal_obs_mea_cov(runs_data);

ln_data = reshape(ln_data, size(ln_data,1)*size(ln_data,2), size(ln_data,3));
ln_data2.values     = ln_data;
ln_data2.model_name = model_name_large;
iv_cov              = cal_obs_iv_cov(ln_data2);
obs_cov = obs_mea_cov + iv_cov; 

[nSpace, nTime, nModel_his] = size(mod_his);
mod_his_stacked = reshape(mod_his, [nSpace*nTime, nModel_his]);  
mo_best  = mean(mod_his_stacked , 2, 'omitnan'); 
M_prior  = cal_mo_cov(mod_his_stacked , scale_m);
Tall =  length(mo_best);

%% run KCC 

Nres=99;
post_all = kriging(obs_best, obs_cov, mo_best, M_prior, (1:Tall)', Nres); 

save('F:\data_code\output\tas_global_GSATconstraint_1973_2024.mat', 'post_all','obs_ref');

