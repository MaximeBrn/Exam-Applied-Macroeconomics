function [var_mean,var_CI_lower,var_CI_upper] = get_CI(IRFs)
    % Note:
    % quantiles(matrix,p,2) -> sum over the raws
    % So, i-th raw quantiles(matrix,p,2) = p-quantile response value at the i-th quarter 
    var_mean=quantile(IRFs,0.5,2); 
    var_CI_lower=quantile(IRFs,0.05,2); 
    var_CI_upper=quantile(IRFs,0.95,2);
end
