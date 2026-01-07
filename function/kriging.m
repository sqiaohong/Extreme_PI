%% ================== Kriging ==================

function post_series = kriging(obs_value, obs_cov, mod_value, mod_cov, time_vec, Nres)
% obs_value : y [T_his x 1]
% obs_cov   : R  [T_his x T_his]
% mod_value : m[T_all x 1]
% mod_cov   : M  [T_all x T_all]
% time_vec  : [T_all x 1]（
% Nres      : default 99
    if nargin < 6 || isempty(Nres), Nres = 99; end

    This = numel(obs_value);
    Tall = numel(mod_value);
    H    = H_matrix(This, Tall);

    S    = H*mod_cov*H' + obs_cov; 
    Sinv = pinv(S); 
    K    = mod_cov * H' * Sinv; 
 
    x_post  = mod_value + K * (obs_value - H*mod_value);  
    cov_post = mod_cov - K*H*mod_cov; 
 
   
    Msqrt = matrix_sqrt_eigen(cov_post);  
    ns    = Tall; 
    rng(42, 'twister'); 
    epsil = zeros(ns, Nres+1);
    epsil(:,2:end) = randn(ns, Nres);

    post_series = repmat(x_post, 1, Nres+1) + Msqrt * epsil;  

end

