function [ phi ] = compute_phi(  )
    %FUNCTION COMPUTE_PHI computes integral phi required for Fisher information
    %   For large x, integral is difficult to compute numerically 
    %   Thus we use the fact that phi is approximately 1 for large x
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014 
    
     % function to integrate with respect to y
     f =  @(y, z) 1/z*y^2 * log(besseli(0, y*z)) * exp(-0.5*(z^2+y^2)) ...
         *((y^2 - 3)*besseli(1, z*y) ...
         - 0.5*y*z*(besseli(0, z*y) + besseli(2, z*y))); 

     % upper limit of integration 
     % must be chosen small enough that besseli does not return infinity
     % but large enough that integral from ymax to infinity is negligible 
     ymax = 50; 
     
     % points z at which to integrate 
     zvals = linspace(0.00000001, 20, 500); 
     
     % set number of integrals to compute
     % must be chosen large enough that phi(zvals(index_max)) has almost 
     % converged to 1 but small enough that integral from ymax to infinity
     % is negligible     
     index_max = 350; 
     
     % compute phi
     phi_vals = zeros(length(zvals), 1); 
     for i = 1:length(zvals)
         z = zvals(i); 
         h = @(y) f(y , z); 
         if i <= index_max
           % compute phi by computing integral from 0 to ymax 
           phi_vals(i) = integral(h, 0, ymax, 'ArrayValued', true) ...
               - zvals(i)^2; 
         else
           % for values of z too large, just set phi(z) = 1
           phi_vals(i) = 1; 
         end
     end
     
     % plot phi_vals to check that output is correct 
     % plot(zvals, phi_vals)
     % xlabel('SNR')
     % ylabel('\phi')

     % create function that interpolates between comput
     phi = griddedInterpolant(zvals, phi_vals); 
         
end

