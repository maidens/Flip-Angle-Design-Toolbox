% set up tests
function tests = least_squares_objective_test
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  path(path, '../')
end

function teardownOnce(testCase)
  rmpath('../')
end

% set tolerance for tests 
function eps = epsilon()
  eps = 1e-5;
end

function testEqualOutputs(testCase)

    % generate random 2D model
    A = randn(2, 2) - eye(2); 
    B = randn(2, 1); 
    C = eye(2); 
    D = [0; 0]; 
    
    % discretize the model using zero-order hold 
    TR = 1.5; 
    sys = ss(A, B, C, D); 
    sysd = c2d(sys, TR); 
    
    % set input, number of samples, initial state
    N = 100; 
    u_fun = @(t) sin(t/10); 
    x0 = [1; 0];
    times = TR*(1:N); 
    for t=1:N
        u(:, t) = u_fun(times(t)); 
    end
    
    % all zero flip angles 
    thetas = zeros(2, N); 
    
    % compute output of trajectories function 
    [y, x] = trajectories(thetas, sysd.a, sysd.b, sysd.c, sys.d, u, x0, N); 
    
    % compute true value using MATLAB's lsim function with zero-order hold
    %    discretization 
    y_lsim = lsim(sys, u_fun(times), times, x0, 'zoh');  
    
    % check that xtilde is zero, and x is the same as y, u is the same as u_fun 
    verifyEqual(testCase, x', y_lsim, 'AbsTol', epsilon())
    verifyEqual(testCase, y', zeros(size(y')), 'AbsTol', epsilon())
    verifyEqual(testCase, u, u_fun(times), 'AbsTol', epsilon())

end
