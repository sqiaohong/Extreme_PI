function [x_vals, P_sorted] = ecdf_plot_inverse(data, P_targets)

data = data(~isnan(data)); 
    data_sorted = sort(data);  
    n = length(data_sorted);
    P_sorted = ((1:n) - 0.5) / n;  
    x_vals = interp1(P_sorted, data_sorted, P_targets, 'linear', 'extrap');
    
end
