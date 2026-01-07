function [Sigma_shrink, shrinkage] = ledoit_wolf_cov(X)
% Ledoit–Wolf shrinkage covariance estimator

    [n, p] = size(X);
    X = X - mean(X, 1);         
    S = (X' * X) / n;                

    mu = trace(S) / p;               
    F = mu * eye(p);                 
    
    X2 = X.^2;
    phiMat = (X2' * X2) / n - 2 * (S .* S) + S.^2;
    phi = sum(phiMat(:));
    rho = sum((S(:) - F(:)).^2);
    kappa = phi / rho;
    shrinkage = max(0, min(1, kappa / n));

    Sigma_shrink = shrinkage * F + (1 - shrinkage) * S;
end
