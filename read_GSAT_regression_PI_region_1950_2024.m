clear
clc

load('F:\data_code\output\tas_global_GSATconstraint_1950_2024.mat')
load('F:\data_code\input\ERA5_global_ecdf_1950_2025.mat')

year=1950:2024;
idd1=find(year==1981);
idd2=find(year==2010);
base=mean(series_ERA5(idd1:idd2));
ERA5=series_ERA5(1:end-1)-base;
tas = obs_ref-mean(obs_ref(idd1:idd2));
post_all=post_all-mean(post_all(idd1:idd2,:),1);
ALL=post_all(:,1);
ALL_train = ALL(1:end-1);
NAT=tas-ALL;
NAT_train=NAT(1:end-1);
X = [ALL_train,NAT_train];

pr_train  = ERA5(1:end-1);
y = pr_train;
mdl = fitlm(X, pr_train);
beta = mdl.Coefficients.Estimate;
CovB = mdl.CoefficientCovariance;
ALL_2024=ALL(end);
NAT_2024=NAT(end);
x_2024 = [1, ALL_2024,NAT_2024] ;
yhat_2024_global = x_2024 * beta;
yhat_ALL_global = beta(2) *ALL + beta(1)/2;
yhat_NAT_global = beta(3) *NAT + beta(1)/2;


%% region result
load('F:\data_code\input\ERA5_region_ecdf_1950_2025.mat')
series_ERA5_region=series_ERA5_region(1:end-1,:);
base=mean(series_ERA5_region(idd1:idd2,:),1);
series_ERA5_region=series_ERA5_region-base;

for m=1:44
    pr_train  = series_ERA5_region(1:end-1,m);
    y = pr_train;
    mdl = fitlm(X, pr_train);
    beta = mdl.Coefficients.Estimate;
    CovB = mdl.CoefficientCovariance;
    ALL_2024=ALL(end);
    NAT_2024=NAT(end);
    x_2024 = [1, ALL_2024,NAT_2024] ;
    yhat_2024_region(m) = x_2024 * beta;

    yhat_ALL_region(:,m)  = beta(2) *ALL + beta(1)/2;
    yhat_NAT_region(:,m)  = beta(3) *NAT + beta(1)/2;

end

save('F:\data_code\output\GSAT_regression_PI_region_1950_2024.mat','yhat_NAT_region', 'yhat_ALL_region', 'yhat_2024_region', 'yhat_NAT_global', 'yhat_ALL_global', 'yhat_2024_global')

