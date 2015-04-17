% set up tests
function tests = compute_phi_test
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  path(path, '../')
  testCase.TestData.phi = compute_phi(); 
end

function teardownOnce(testCase)
  clear testCase.TestData.phi
  rmpath('../')
end

% set tolerance for tests 
function eps = epsilon()
  eps = 1e-5;
end

% test that phi is zero at zero
function testLowerLimit(testCase)
    verifyEqual(testCase, testCase.TestData.phi(0), 0, 'AbsTol', epsilon())
end

% test that phi goes to 1 as SNR goes to infinity 
function testUpperLimit(testCase)
    bigNumber = 1e10; 
    verifyEqual(testCase, testCase.TestData.phi(bigNumber), 1, 'AbsTol', epsilon())
end

% compare results of interpolation in phi with direct computation 
function testInterpolation(testCase)
    % values at which to compute integral 
    svals = logspace(-1, 1, 10); 
    
    % function to integrate 
    function I = integrand(s, x) 
        I = x^3*besseli(1,s*x,1)^2/besseli(0,s*x,1)*exp(-(s^2 -2*s*x+ x^2)/2); 
        % This is equal to x^3*besseli(1,s*x)^2/besseli(0,s*x)*exp(-(s^2 + x^2)/2)
        %   but scales better for large s*x 
        %   see MATLAB documentation for besseli(nu, Z, 1)
    end

    % compute integral wrt x at points s in svals 
    test_phi_vals = zeros(length(svals), 1); 
    for i = 1:length(svals)
         h = @(x) integrand(svals(i), x); 
         % compute phi by computing integral from 0 to Inf 
         test_phi_vals(i) = integral(h, 0, Inf, 'ArrayValued', true) - svals(i)^2; 
    end
    
    % compute values of phi to compare 
    phi_vals = zeros(length(svals), 1); 
    for i = 1:length(svals)
         phi_vals(i) = testCase.TestData.phi(svals(i)); 
    end
    
    % check that this agrees with phi 
    verifyEqual(testCase, phi_vals, test_phi_vals, 'AbsTol', epsilon())
end
    