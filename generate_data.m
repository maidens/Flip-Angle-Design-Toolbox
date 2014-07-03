function [y, y_true] = generate_data(model, thetas)
    %FUNCTION GENERATE_DATA computes simulated trajectories of the model 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014     
    
    
    % evaluate symbolic variables at values provided 
    Ad_nom = double(subs(model.Ad_sym, [model.parameters_of_interest, ...
        model.nuisance_parameters, model.known_parameters], ...
        [model.parameters_of_interest_nominal_values, ...
        model.nuisance_parameters_nominal_values, ...
        model.known_parameter_values]));
    Bd_nom = double(subs(model.Bd_sym, [model.parameters_of_interest, ...
        model.nuisance_parameters, model.known_parameters], ...
        [model.parameters_of_interest_nominal_values, ...
        model.nuisance_parameters_nominal_values, ...
        model.known_parameter_values]));
    u_fun = matlabFunction(subs(model.u,[model.parameters_of_interest, ...
        model.nuisance_parameters, model.known_parameters], ...
        [model.parameters_of_interest_nominal_values, ...
        model.nuisance_parameters_nominal_values, ...
        model.known_parameter_values]));
   
    % compute "true" trajectories of model 
    y_true = trajectories(thetas, Ad_nom, Bd_nom, u_fun, model.TR, model.N);
    
    if strcmp(model.noise_type, 'None')
        
        % add no noise to data 
        y = y_true;
        
    elseif strcmp(model.noise_type, 'Rician')
        
        % add Rician noise to data 
        sigma_vector = sqrt(model.noise_parameters); 
        r1 = randn(size(y_true));
        r2 = randn(size(y_true));
        y = sqrt((y_true + diag(sigma_vector)*r1).^2 ...
            + (diag(sigma_vector)*r2).^2);
        
    else
        error(strcat('Noise type ', model.noise_type, 'is invalid'))
    end
    
end
   



