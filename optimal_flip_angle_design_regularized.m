function [thetas_opt, objective_value, q_opt] = optimal_flip_angle_design_regularized(model, ...
    design_criterion, initial_q_value, lambda, options)
    %FUNCTION OPTIMAL_FLIP_ANGLE_DESIGN performs numerical optimization to generate optimal flip angle scheme
    %   Choice of design criteria: 
    %       * 'totalSNR'  (maximize sum of state and input variables over all time) 
    %       * 'D-optimal' (maximize determinant of information matrix)
    %       * 'E-optimal' (maximize smalles eigenvalues of information matrix)
    %       * 'T-optimal' (maximize trace of information matrix)
    %       * 'A-optimal' (minimize trace of inverse information matrix)
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
    % bounds to keep flip angles between 0 and 90 degrees 
    lb = zeros(size(initial_q_value)); 
    ub = [pi/2*ones(size(initial_q_value))]; 
    % ub = pi/2*ones(size(initial_q_value)); 
    

    display('===== Computing optimal flip angles =====')
       
    % Optimal flip angle design for total SNR design criterion 
    if strcmp(design_criterion,'totalSNR')
        % define objective function for total SNR 
        obj = @(q) -sum(sum(trajectories(model.flip_angle_input_matrix*q, model.Ad_nom, ...
            model.Bd_nom, model.C_nom, model.D_nom, model.u_fun(model.TR*0:model.N), model.x0_nom, model.N)));

        % initialize optimization problem 
        if nargin < 4
            options = optimset('MaxFunEvals', 50000, 'MaxIter', 1000, ...
                'Display', 'iter'); 
        end
        %initial_thetas_value = ones(model.N, 3); 
        
        % perform optimization 
        [q_opt, objective_value] = fmincon(obj, initial_q_value, ...
                             [], [], [], [], lb, ub, [], options);
        thetas_opt = model.flip_angle_input_matrix * q_opt; 
        
    % Optimal flip angle design for Fisher information design criteria
    elseif strcmp(design_criterion,'D-optimal') || ...
            strcmp(design_criterion,'E-optimal') || ...
            strcmp(design_criterion,'A-optimal') || ...
            strcmp(design_criterion,'T-optimal')
        
        if strcmp(model.noise_type,'Rician')
            
            % compute function phi 
            phi = compute_phi(); 
            
            % check that model structure is identifiable 
            information = det(fisher_information(model.flip_angle_input_matrix*initial_q_value, model, phi)); 
            if information == 0 || isnan(information)
                error('Fisher information matrix is singular. Try ensuring that model structure (A, B, u) = (%s, %s, %s) is identifiable for parameter vector p = [%s %s] or choose a different initial value for the flip angles.', char(model.A), char(model.B), char(model.u), char(model.parameters_of_interest), char(model.nuisance_parameters)) 
            end

            % define objective function 
            if strcmp(design_criterion,'D-optimal')
                obj_coarse = @(q) ...
                    -log(abs(det(fisher_information(model.flip_angle_input_matrix*q, model, phi)))) + lambda*norm(reshape(q(:, 2:end) - q(:, 1:end-1), 2*model.N-2, 1), 2)   ; 
            elseif strcmp(design_criterion,'E-optimal')
                obj_coarse = @(q) ...
                    -min(eig(fisher_information(model.flip_angle_input_matrix*q, model, phi)));
            elseif strcmp(design_criterion,'T-optimal')
                obj_coarse = @(q) ...
                    -trace(fisher_information(model.flip_angle_input_matrix*q, model, phi));
            elseif strcmp(design_criterion,'A-optimal')
                obj_coarse = @(q) ...
                    trace(inv(fisher_information(model.flip_angle_input_matrix*q, model, phi)));
            else
                error('this should not be reachable -- possibly error with design_criterion handling')
            end
            
            % set options for optimization problem 
            if nargin < 4
                options = optimset('MaxFunEvals', 5000, 'MaxIter', 200, ...
                    'Display', 'iter'); 
            end
            
            % solve optimization problem 
            [q_opt, objective_value] = fmincon(obj_coarse, initial_q_value, ...
                 [], [], [], [], lb, ub, [], options);
            thetas_opt = model.flip_angle_input_matrix * q_opt;            
        else
            
            error(['Design criterion "', design_criterion, ...
                '" is not valid for noise of type "', ...
                model.noise_type, '"'])
        end
                
    else
        error(['Design criterion "', design_criterion, '" is invalid'])
    end
    
end




