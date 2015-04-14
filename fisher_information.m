function [I_Schur, I_full] = fisher_information(thetas, model, phi)
    %FUNCTION FISHER_INFORMATION Computes information matrix for a particular flip angle scheme 
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
    u = zeros(model.ni, model.N); 
    for t=0:model.N-1
        u(:,t+1)  = model.u_fun(model.TR*t); 
    end
    [y, x] = trajectories(thetas, model.Ad_nom, model.Bd_nom, model.C_nom, ...
        model.D_nom, u, model.x0_nom, model.N); 
    
    % compute sensitivities of state with respect to parameters
    p = [model.parameters_of_interest, model.nuisance_parameters];   
    DuDp = zeros(model.ni, model.N, length(p)); 
    DxDp = zeros(model.n, model.N, length(p));
    DyDp = zeros(model.no, model.N, length(p)); 
    
    % loop over parameters 
    for i=1:length(p)
        
        % unpack system sensitivites 
        DAi = model.sensitivity_Ad(:, :, i);  
        DBi = model.sensitivity_Bd(:, :, i);
        DCi = model.sensitivity_C(:, :, i); 
        DDi = model.sensitivity_D(:, :, i);
        DuDp(:, :, i) = model.sensitivity_u(:, :, i);
        
        % compute state trajectory sensitivities 
        DxDp(:, 1, i) = model.sensitivity_x0(:, i);
        for t=1:model.N-1
            DxDp(:, t+1, i) = ...
                DAi * diag(cos(thetas(:, t))) * x(:, t) ...
                + model.Ad_nom * diag(cos(thetas(:, t))) * DxDp(:, t, i) ...
                + DBi * u(:, t) ...
                + model.Bd_nom * DuDp(:, t, i); 
        end
        
        % compute output trajectory sensitivities 
        DyDp(:, 1, i) = DCi*x(:, 1) + model.C_nom*DxDp(:, 1, i); 
        for t=1:model.N 
            DyDp(:, t, i) = ...
                DCi * diag(sin(thetas(:, t))) * x(:, t) ...
                + model.C_nom * diag(sin(thetas(:, t))) * DxDp(:, t, i) ...
                + DDi * u(:, t) ...
                + model.D_nom * DuDp(:, t, i);  
        end

    end
    
    % compute integral phi along trajectories 
    phi_of_x = [];
    for k = 1:model.no
        phi_of_x = [phi_of_x; 
                    phi(y(k, :) ...
                    /sqrt(model.noise_parameters(k)))];
    end
    
    % compute information matrix 
    I_full = zeros(length(p)); 
    for i=1:length(p)
        for j=1:i
            I_full(i, j) = sum(sum( ...
                diag(model.noise_parameters.^(-1)) ...
                * (DyDp(:, :, i) .* DyDp(:, :, j) .* phi_of_x) ...
                )); 
            I_full(j, i) = I_full(i, j); 
        end
    end

    % compute Schur complement corresponding to parameters of interest 
    len = length(model.parameters_of_interest); 
    I_Schur = I_full(1:len, 1:len) ...
       - I_full(1:len, len+1:end)*(I_full(len+1:end, len+1:end)\I_full(len+1:end, 1:len)); 

end

