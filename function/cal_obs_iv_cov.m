function lw_cov = cal_obs_iv_cov(ln_data)
% cal_obs_iv_cov - Compute internal variability covariance using Ledoit–Wolf shrinkage

    % -------------------------------
    if istable(ln_data)
        values = ln_data.values;
        model_name = ln_data.model_name;
    elseif isstruct(ln_data)
        values = ln_data.values;
        model_name = ln_data.model_name;
    else
        error('ln_data 必须是 struct 或 table，包含 fields: values, model_name');
    end


    ln_name = unique(model_name, 'stable');
    nModel  = numel(ln_name);

    lw_each = cell(nModel, 1);

    % -------------------------------
    % Ledoit–Wolf 
    % -------------------------------
    for im = 1:nModel
   
        idx = strcmp(model_name, ln_name{im});
        mod_data = values(:, idx);  
        

        lw_each{im} = ledoit_wolf_cov(mod_data');
    end


    lw_cov = mean(cat(3, lw_each{:}), 3);

end


