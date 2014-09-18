function [parameters_of_interest_est, nuisance_parameters_est] ...
    = simultaneous_ML_estimation(model, y1, y2, thetas1, thetas2)
    % FUNCTION PARAMETER_ESTIMATION computes an estimate of model parameters from data y
    %   Choice of goodness_of_fit_criterion: 
    %      * 'least-squares'
    %      * 'maximum-likelihood'
    %
    % Flip Angle Design Toolbox 
    % John Maidens (maidens@eecs.berkeley.edu) 
    % June 2014 
    
    display('===== Estimating parameter values =====') 
   
    % perform maximum-likelihood fit 
    options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000); 
    obj = @(p) negative_log_likelihood_rician(p, y1, ...
        thetas1*model.flip_angle_input_matrix', model) ...
     + negative_log_likelihood_rician(p, y2, ...
        thetas2*model.flip_angle_input_matrix', model); 

    % initial point 
    p_init =  [model.parameters_of_interest_nominal_values, ...
        model.nuisance_parameters_nominal_values]; 

    p_opt = fminunc(obj, p_init, options);

    % unpack estimates 
    l = length(model.parameters_of_interest_nominal_values); 
    parameters_of_interest_est = p_opt(1:l); 
    nuisance_parameters_est = p_opt(l+1:end); 

end
