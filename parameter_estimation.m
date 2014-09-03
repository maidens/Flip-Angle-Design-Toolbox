function [parameters_of_interest_est, nuisance_parameters_est] ...
    = parameter_estimation(y, model, goodness_of_fit_criterion, thetas)
    % FUNCTION PARAMETER_ESTIMATION computes an estimate of model parameters from data y
    %   Choice of goodness_of_fit_criterion: 
    %      * 'least-squares'
    %      * 'maximum-likelihood'
    %
    % Flip Angle Design Toolbox 
    % John Maidens (maidens@eecs.berkeley.edu) 
    % June 2014 
    
    display('===== Estimating parameter values =====')
    
    if strcmp(goodness_of_fit_criterion, 'least-squares')
        
        % perform least-squares fit 
        options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000); 
        obj = @(p) least_squares_objective(p, y, thetas*model.flip_angle_input_matrix', model); 
        p_opt = lsqnonlin(obj, ...
            [model.parameters_of_interest_nominal_values, ...
            model.nuisance_parameters_nominal_values], [], [], options);
        
        % unpack estimates 
        l = length(model.parameters_of_interest_nominal_values); 
        parameters_of_interest_est = p_opt(1:l); 
        nuisance_parameters_est = p_opt(l+1:end); 
        
        
    elseif strcmp(goodness_of_fit_criterion, 'maximum-likelihood')
        if strcmp(model.noise_type, 'Rician')
            
            % perform maximum-likelihood fit 
            options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000); 
            %obj = @(p) negative_log_likelihood_rician(p, y, thetas*model.flip_angle_input_matrix', model); 
            obj = @(p) negative_log_likelihood_rician(p, y, thetas*model.flip_angle_input_matrix', model); 
            p_opt = fminunc(obj, ...
                [model.parameters_of_interest_nominal_values, ...
                model.nuisance_parameters_nominal_values], options);

            % unpack estimates 
            l = length(model.parameters_of_interest_nominal_values); 
            parameters_of_interest_est = p_opt(1:l); 
            nuisance_parameters_est = p_opt(l+1:end); 
        
        else
            error(['Cannot do maximum likelihood fit for noise of type "', ...
                model.noise_type ,'"'])
        end
        
        
        
    else
        error(['Goodness of fit criterion "' , goodness_of_fit_criterion , ...
            '" is invalid'])
    end
    
    
end
