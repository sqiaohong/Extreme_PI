function [post_all, prior, mo_best, M_prior, obs_best, obs_cov] = ...
    run_kcc_constraint_1d(obs_ref, obs_200runs_ref, ln_GSAT_data_ref, mod_ref, ...
    mo_da_fu, model_name_large, model_GSAT_name, model_PI, ...
    scale_m, Nres)

%   obs_ref            : (Tobs × 1)  
%   obs_200runs_ref    : (Tobs × Nrun) 
%   ln_GSAT_data_ref   : (Tobs × Niv)   
%   mod_ref            : (Tprior × Nmodel_his) 
%   mo_da_fu           : (Tfuture × Nmodel_PI) 
%   model_name_large   : Niv×1 cell
%   model_fut        
%   model_PI          
%   scale_m          
%   Nres      
%

% -------- 1. select obs and model data for constraint --------
obs_best          = obs_ref(:);        
runs_data(:,1,:)  = obs_200runs_ref;  
ln_data(:,1,:)    = ln_GSAT_data_ref; 
mod_his(:,1,:)    = mod_ref;          

% -------- 2. obs convariate R = measure error + internal variability --------
runs_data = reshape(runs_data, size(runs_data,1)*size(runs_data,2), size(runs_data,3));
obs_mea_cov = cal_obs_mea_cov(runs_data);

ln_data = reshape(ln_data, size(ln_data,1)*size(ln_data,2), size(ln_data,3));
ln_data2.values     = ln_data;
ln_data2.model_name = model_name_large;
iv_cov              = cal_obs_iv_cov(ln_data2);

obs_cov = obs_mea_cov + iv_cov;   

% -------- 3. model prior --------
[nSpace, nTime, nModel_his] = size(mod_his);
mod_his_stacked = reshape(mod_his, [nSpace*nTime, nModel_his]); 

% 
[common_models, ia, ib] = intersect(model_GSAT_name, model_PI, 'stable');
mod_his_stacked1 = mod_his_stacked(:, ia); 
mo_da_fu_sel     = mo_da_fu(:, ib);        

%
mo_his_da_or_fu = cat(1, mod_his_stacked1, mo_da_fu_sel);  
nan_cols       = any(isnan(mo_his_da_or_fu), 1);
mo_his_da_or_fu = mo_his_da_or_fu(:, ~nan_cols);

% prior MME and convariate
mo_best = mean(mo_his_da_or_fu, 2, 'omitnan');    
M_prior = cal_mo_cov(mo_his_da_or_fu, scale_m);    

% -------- 4. H matrix --------
Thobs = length(obs_best(:));
Tall  = length(mo_best);
H     = H_matrix(Thobs, Tall); 

% -------- 5. Kriging--------
post_all_full = kriging(obs_best, obs_cov, mo_best, M_prior, (1:Tall)', Nres);

Tfu     = size(mo_da_fu_sel, 1);
post_all = post_all_full(end-Tfu+1:end, :); 
prior    = mo_da_fu_sel(:, ~nan_cols);      

end
