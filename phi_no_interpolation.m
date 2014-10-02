function [ phi ] = phi_no_interpolation(  )
    %FUNCTION COMPUTE_PHI computes integral phi required for Fisher information
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
     % function to integrate with respect to y
     f = @(y, z) y^3*besseli(1,z*y,1)^2/besseli(0,z*y,1)*exp(-(z^2 -2*z*y+ y^2)/2);
     % This is equal to x^3*besseli(1,s*x)^2/besseli(0,s*x)*exp(-(s^2 + x^2)/2)
        %   but scales better for large s*x 
        %   see MATLAB documentation for besseli(nu, Z, 1)

     phi = @(z) integral(@(y)f(y, z), 0, Inf, 'ArrayValued', true) - z.^2; 
end

