function [ phi ] = compute_phi(  )
    %FUNCTION COMPUTE_PHI computes integral phi required for Fisher information
    %   Since phi will be called many times, we compute its value at points
    %   on a logarithmically spaced grid then use interpolation to evaluate
    %   it when needed. This provides significant speedup in the
    %   optimization. For large x, integral is difficult to compute 
    %   numerically so we use the fact that phi is approximately 1 for 
    %   large x
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
     % function to integrate with respect to y
     f = @(y, z) y^3*besseli(1,z*y,1)^2/besseli(0,z*y,1)*exp(-(z^2 -2*z*y+ y^2)/2);
     % This is equal to x^3*besseli(1,s*x)^2/besseli(0,s*x)*exp(-(s^2 + x^2)/2)
        %   but scales better for large s*x 
        %   see MATLAB documentation for besseli(nu, Z, 1)
        
     
     % points z at which to integrate 
     % zvals = linspace(0.00000001, 20, 500); 
     num_z = 1000; 
     zvals = [0 logspace(-2, 2, num_z)]; 
     
     % for i > index_max we set phi_vals equal to 1 to ensure that
     %      lim_z->Inf phi(z) = 1 
     index_max = num_z-10; 
     
     % compute phi
     phi_vals = zeros(length(zvals), 1); 
     for i = 1:length(zvals)
         z = zvals(i); 
         h = @(y) f(y , z); 
         if i <= index_max
           % compute phi by computing integral from 0 to ymax 
           phi_vals(i) = integral(h, 0, Inf, 'ArrayValued', true) ...
               - zvals(i)^2; 
         else
           % for large, just set phi(z) = 1
           phi_vals(i) = 1; 
         end
     end

     % create function that interpolates between comput
     phi = griddedInterpolant(zvals, phi_vals); 
         
end

