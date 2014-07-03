function obj = least_squares_objective(p, y, thetas, model)
    %FUNCTION LEAST_SQUARES_OBJECTIVE computes value of objective funtion for least-squares parameter fitting  
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014     
    
    % set value of model parameters to values passed in argument p
    l = length(model.parameters_of_interest_nominal_values);
    model.parameters_of_interest_nominal_values = p(1:l);
    model.nuisance_parameters_nominal_values = p(l+1:end); 
    
    % compute model trajectories with given parameter values
    [~, y_true] = generate_data(model, thetas); 
    
    % return difference between measured data and computed trajectories 
    obj = y - y_true;
    
end