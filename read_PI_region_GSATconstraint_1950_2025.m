
%%
clear
clc

addpath('F:\data_code\input\', 'F:\data_code\output\', 'F:\data_code\function\');
start_ref=1950;
end_ref=2024;

start_year=1950;
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


%% Smoothed PI from 1950-2025
load('Rx1day_ecdf_25dg_lon180_landmean_noAntarctica_1950_2025.mat')
load('model_region_ecdf_1950_2025.mat')
[model_PI, ia,ib] = intersect(model_PI, PI_model_name, 'stable');
region_PI_model=squeeze(series_model_region(:,:,ib));
year_PI=start_year:end_year+1;
for m=1:size(region_PI_model,2)
    for n=1:size(region_PI_model,3)
        data=region_PI_model(:,m,n);
        smooth1 = ncs_random_smooth(year_PI, data, 'NRepeat', 10, 'SpanYears', [5 15], 'ReturnAll', true);
        PI_smooth(:,m,n)=smooth1.y_smooth;
    end
end
iy=find(year_PI>=start_ref & year_PI<=end_ref);
mo_da_fu = PI_smooth-mean(PI_smooth(iy,:,:),1); 


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

%% run KCC
scale_m=1;
Nres=99;
for m=1:size(region_PI_model,2)
    
    mo_da_fu_region=squeeze(mo_da_fu(:,m,:));

    %% KCC
    [post_all_region, prior_region, mo_best_region, M_prior_region, obs_best_region, obs_cov_region] = ...
        run_kcc_constraint_1d(obs_ref, obs_200runs_ref, ln_GSAT_data_ref, mod_ref, ...
        mo_da_fu_region, model_name_large, model_GSAT_name, model_PI, ...
        scale_m, Nres);

    %% save data
    post_all(:,:,m) = post_all_region;
    prior(:,:,m)=prior_region;

end

save('F:\data_code\output\PI_region_GSATconstraint_1950_2025.mat', 'post_all', 'prior');

