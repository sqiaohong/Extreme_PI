clear
clc

load('F:\data_code\output\tas_global_GSATconstraint_1973_2024.mat')
load('F:\data_code\input\ERA5_global_ecdf_1973_2025.mat')

year=1973:2024;
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




save('F:\data_code\output\GSAT_regression_PI_region_1973_2024.mat','yhat_NAT_global', 'yhat_ALL_global', 'yhat_2024_global')

