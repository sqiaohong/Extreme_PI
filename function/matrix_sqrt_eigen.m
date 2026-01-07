function Msqrt = matrix_sqrt_eigen(M)
    [V, D] = eig((M+M')/2);              
    ev     = max(real(diag(D)), 0);
    Msqrt  = real(V * diag(sqrt(ev)) / V);
end
