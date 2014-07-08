function I_Schur = fisher_information(thetas, model, phi)
    %FUNCTION FISHER_INFORMATION Computes information matrix for a particular flip angle 
    %   Returns Schur complement that corresponds to eliminating parameters
    %   stored in model.nuisance_parameters
    %   
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
        
    % compute discretized model (if necessary) 
    if ~model.discretized 
        model = discretize(model); 
    end
    
    % compute system matrix and input sensitivities (if necessary) 
    if ~model.sensitivities_computed
        model = sensitivities(model); 
    end
    
    % get input and state trajectories 
    y = trajectories(thetas, model.Ad_nom, model.Bd_nom, model.u_fun, ...
        model.x0_nom, model.TR, model.N);
    x = y(1:model.n, :);
    u_val = y(model.n+1:model.n+model.m, :);
    
    % compute sensitivities of state with respect to parameters
    p = [model.parameters_of_interest, model.nuisance_parameters];   
    DuDp = zeros(model.N, length(p)); 
    DxDp = zeros(model.n, model.N, length(p)); 
    
    for i=1:length(p)
        DAi = model.sensitivity_Ad(:, :, i);  
        DBi = model.sensitivity_Bd(:, :, i); 
        DxDp(:, 1, i) = model.sensitivity_x0(:, i); 
        for t=1:model.N-1
            DuDp(t, i) = model.sensitivity_u(t, i);
            if (t < model.N)
                DxDp(:, t+1, i) = ...
                    DAi * diag(cos(thetas(t, 1:model.n))) * x(:, t) ...
                    + model.Ad_nom * diag(cos(thetas(t, 1:model.n))) * DxDp(:, t, i) ...
                    + DBi * u_val(t) ...
                    + model.Bd_nom * DuDp(t, i); 
            end
        end
    end
    
    % compute integral phi along trajectories 
    phi_of_x = [];
    for i = 1:(model.n + model.m)
        phi_of_x = [phi_of_x; 
                    phi(sin(thetas(:, i))'.*y(i,:) ...
                    /sqrt(model.noise_parameters(i)))];
    end
    
    % compute information matrix 
    I = zeros(length(p)); 
    for i=1:length(p)
        for j=i:length(p)
            I(i, j) = sum(sum( sin(thetas)'.^2 ...
                .*(diag(model.noise_parameters.^-1) ...
                * ([DxDp(:, :, i); DuDp(:, i)'] .* phi_of_x ...
                .* [DxDp(:, :, j); DuDp(:, j)'])))); 
            I(j, i) = I(i, j); 
        end
    end

    % compute Schur complement corresponding to parameters of interest 
    len = length(model.parameters_of_interest); 
    I_Schur = I(1:len, 1:len) ...
       - I(1:len, len+1:end)*(I(len+1:end, len+1:end)\I(len+1:end, 1:len)); 

end

