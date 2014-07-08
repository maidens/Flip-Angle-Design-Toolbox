function thetas_opt = optimal_flip_angle_design(model, ...
    design_criterion, initial_thetas_value, options)
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
       
    % Optimal flip angle design for total SNR design criterion 
    if strcmp(design_criterion,'totalSNR')
        % define objective function for total SNR 
        obj = @(thetas) -sum(sum(trajectories(thetas, model.Ad_nom, ...
            model.Bd_nom, model.u_fun, model.x0_nom, model.TR, model.N)));
        
        % initialize optimization problem 
        if nargin < 4
            options = optimset('MaxFunEvals', 50000, 'MaxIter', 1000, ...
                'Display', 'iter'); 
        end
        %initial_thetas_value = ones(model.N, 3); 
        
        % perform optimization 
        thetas_opt = fminunc(obj, initial_thetas_value, options);
          
        
    % Optimal flip angle design for Fisher information design criteria
    elseif strcmp(design_criterion,'D-optimal') || ...
            strcmp(design_criterion,'E-optimal') || ...
            strcmp(design_criterion,'T-optimal')
        
        if strcmp(model.noise_type,'Rician')
            
            % compute function phi 
            phi = compute_phi(); 
            
            % check that model structure is identifiable 
            information = det(fisher_information(initial_thetas_value, model, phi)); 
            if information == 0 || isnan(information)
                error('Fisher information matrix is singular. Try ensuring that model structure (A, B, u) = (%s, %s, %s) is identifiable for parameter vector p = [%s %s] or choose a different initial value for the flip angles.', char(model.A), char(model.B), char(model.u), char(model.parameters_of_interest), char(model.nuisance_parameters)) 
            end
            
            % define objective function 
            if strcmp(design_criterion,'D-optimal')
                obj_coarse = @(thetas) ...
                    -log(abs(det(fisher_information(thetas, model, phi)))); 
            elseif strcmp(design_criterion,'E-optimal')
                obj_coarse = @(thetas) ...
                    -min(eig(fisher_information(thetas, model, phi)));
            elseif strcmp(design_criterion,'T-optimal')
                obj_coarse = @(thetas) ...
                    -trace(fisher_information(thetas, model, phi));
            else
                error('this should not be reachable -- possibly error with design_criterion handling')
            end
            
            % set options for optimization problem 
            if nargin < 4
                options = optimset('MaxFunEvals', 5000, 'MaxIter', 200, ...
                    'Display', 'iter'); 
            end
            
            % solve optimization problem 
            thetas_opt = fminunc(obj_coarse, initial_thetas_value, ...
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




