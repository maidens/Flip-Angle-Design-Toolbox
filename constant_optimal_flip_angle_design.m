function thetas_opt = constant_optimal_flip_angle_design(model, ...
    design_criterion, initial_theta_value, options)
    %FUNCTION OPTIMAL_FLIP_ANGLE_DESIGN performs numerical optimization to generate optimal flip angle scheme
    %   Choice of design criteria: 
    %       * 'totalSNR'  (maximize sum of state and input variables over all time) 
    %       * 'D-optimal' (maximize determinant of information matrix)
    %       * 'E-optimal' (maximize smalles eigenvalues of information matrix)
    %       * 'T-optimal' (maximize trace of information matrix) 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 

    display('===== Computing optimal flip angles =====')
       
        
    % Optimal flip angle design for Fisher information design criteria
    if strcmp(design_criterion,'D-optimal') || strcmp(design_criterion,'E-optimal') ||strcmp(design_criterion,'A-optimal') 
           
        if strcmp(model.noise_type,'Rician')
            
            % compute function phi 
            phi = compute_phi(); 
            
            % check that model structure is identifiable 
            information = det(fisher_information(initial_theta_value*ones(model.N, model.n + model.m), model, phi)); 
            if information == 0 || isnan(information)
                error('Fisher information matrix is singular. Try ensuring that model structure (A, B, u) = (%s, %s, %s) is identifiable for parameter vector p = [%s %s] or choose a different initial value for the flip angles.', char(model.A), char(model.B), char(model.u), char(model.parameters_of_interest), char(model.nuisance_parameters)) 
            end
            
            % define objective function 
            if strcmp(design_criterion,'D-optimal')
                obj_coarse = @(theta) ...
                    -log(abs(det(fisher_information(theta*ones(model.N, model.n + model.m), model, phi)))); 
            end
            
            if strcmp(design_criterion,'E-optimal')
                obj_coarse = @(theta) ...
                    -max(eig(fisher_information(theta*ones(model.N, model.n + model.m), model, phi))); 
            end
            
            
            % set options for optimization problem 
            if nargin < 4
                options = optimset('MaxFunEvals', 5000, 'MaxIter', 200, ...
                    'Display', 'iter'); 
            end
            
            % solve optimization problem 
            thetas_opt = fminunc(obj_coarse, initial_theta_value, ...
                options);
                       
        else
            
            error(['Design criterion "', design_criterion, ...
                '" is not valid for noise of type "', ...
                model.noise_type, '"'])
        end
                
    else
        error(['Design criterion "', design_criterion, '" is invalid'])
    end
    
end




