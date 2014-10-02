function [y, y_true] = generate_data(model, thetas)
    %FUNCTION GENERATE_DATA computes simulated trajectories of the model 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014     
    
    
    % unpack system matrices 
    Ad_nom = model.Ad_nom; 
    Bd_nom = model.Bd_nom; 
    C_nom = model.C_nom; 
    D_nom = model.D_nom; 
    u_fun = model.u_fun; 
    u = zeros(model.ni, model.N); 
    for t=1:model.N
        u(:, t) = u_fun(model.TR*(t-1)); 
    end
   
    % compute "true" trajectories of model 
    y_true = trajectories(thetas, Ad_nom, Bd_nom, C_nom, D_nom, u, model.x0_nom, model.N);
    
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
   



