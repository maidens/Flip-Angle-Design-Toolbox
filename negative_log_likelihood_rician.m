function [ l ] = negative_log_likelihood_rician(p, y, thetas, model)
    %FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for compartmental model with Rician noise 
    %   Uses Gaussian approximation to Rician distribution for high SNR 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 

    % compute system trajectory 
    l = length(model.parameters_of_interest_nominal_values);
    model.parameters_of_interest_nominal_values = p(1:l);
    model.nuisance_parameters_nominal_values = p(l+1:end); 
    [~, x_full] = generate_data(model, thetas); 

    % compute negative log likelihood 
    l = 0;
    for t = 1:model.N
        for k = 1:(model.n + model.m)
            % for large arguments, besseli(0, arg) returns inf 
            % thus we use the large argument approximation 
            %   log(besseli(0, x)) ~ x - 0.5*log(2*pi*x) 
            if y(k, t) * x_full(k, t) / model.noise_parameters(k) < 500
                l = l - log(y(k, t)/model.noise_parameters(k)) + (y(k, t)^2 + x_full(k, t)^2) / (2*model.noise_parameters(k)) - log(besseli(0, y(k, t) * x_full(k, t) / model.noise_parameters(k)));
            else
                l = l - log(y(k, t)/model.noise_parameters(k)) + (y(k, t)^2 + x_full(k, t)^2) / (2*model.noise_parameters(k)) - y(k, t) * x_full(k, t) / model.noise_parameters(k) + 0.5*log(2*pi*y(k, t) * x_full(k, t) / model.noise_parameters(k));
            end
        end
    end


end
