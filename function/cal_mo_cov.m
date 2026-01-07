function mo_cov = cal_mo_cov(mo_long_data, scale_m)
%CAL_MO_COV Compute model covariance matrix with scaling
%

    if nargin < 2
        scale_m = 1;
    end

    if isstruct(mo_long_data)
        if isfield(mo_long_data, 'data')
            mo_long_data = mo_long_data.data;
        else
            error('Input struct must contain field "data".');
        end
    end


    len_mod = size(mo_long_data, 2);


    mo_cov = (scale_m^2) * cov(mo_long_data'); 


end
