function [parameters_of_interest_est, nuisance_parameters_est] ...
    = simultaneous_ML_estimation(model, y1, y2, thetas1, thetas2, p_init)

    %
    % Flip Angle Design Toolbox 
    % John Maidens (maidens@eecs.berkeley.edu) 
    % June 2014 
    

        
    % perform maximum-likelihood fit 
    options = optimset('Display','iter','MaxFunEvals', 1000, 'MaxIter', 5000); 
    obj = @(p) negative_log_likelihood_rician(p, y1, ...
        thetas1, model) ...
     + negative_log_likelihood_rician(p, y2, ...
        thetas2, model); 

    % initial point 
    %p_init =  [model.parameters_of_interest_nominal_values, ...
    %    model.nuisance_parameters_nominal_values] 

    p_opt = fminunc(obj, p_init, options);

    % unpack estimates 
    l = length(model.parameters_of_interest_nominal_values); 
    parameters_of_interest_est = p_opt(1:l); 
    nuisance_parameters_est = p_opt(l+1:end); 

end
