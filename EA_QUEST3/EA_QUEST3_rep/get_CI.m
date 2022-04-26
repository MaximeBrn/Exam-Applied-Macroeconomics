function [var_mean,var_CI_lower,var_CI_upper] = get_CI(IRFs)
    var_mean=mean(IRFs,2);
    var_CI_lower=quantile(IRFs,0.05,2); 
    var_CI_upper=quantile(IRFs,0.95,2);
end
