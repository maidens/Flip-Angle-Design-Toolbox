function [y, x] = trajectories(thetas, Ad, Bd, C, D, u, x0, N)
    %FUNCTION TRAJECTORIES computes trajectories of the model 
    %   Not meant to be public, only called from within toolbox functions 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
    n = size(Ad, 1); 
    m = size(C, 1); 
    x = zeros(n, N); 
    y = zeros(m, N); 

    % compute state trajectories 
    x(:, 1) = x0; 
    for t=1:N-1
        x(:, t+1) = Ad * diag(cos(thetas(:, t))) * x(:, t) + Bd * u(:, t);
    end
    % compute output trajectories 
    for t=1:N
        y(:, t) = C * diag(sin(thetas(:, t))) * x(:, t) + D * u(:, t); 
    end
    

end




