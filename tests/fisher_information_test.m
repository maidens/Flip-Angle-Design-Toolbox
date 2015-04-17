% set up tests
function tests = fisher_information_test
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
  eps = 1e-10;
end

% test that fisher_information returns the correct value for a simple model
% with no dynamics for which we can explicitly solve for the Fisher 
% information in terms of phi 
function testStaticModel(testCase)
  
    % set up model with no dynamics
    testCase.TestData.model = linear_exchange_model; 

    syms p

    testCase.TestData.model.A = 0; 
    testCase.TestData.model.B = 0; 
    testCase.TestData.model.C = 1; 
    testCase.TestData.model.D = 0; 
    testCase.TestData.model.u = @(t)(t); 
    testCase.TestData.model.x0 = p; 
    testCase.TestData.model.TR = 1; 
    testCase.TestData.model.N = 1; 

    testCase.TestData.model.parameters_of_interest = p; 
    testCase.TestData.model.parameters_of_interest_nominal_values = 1; 
    testCase.TestData.model.noise_parameters = [1]; 

    testCase.TestData.model = discretize(testCase.TestData.model); 
    testCase.TestData.model = sensitivities(testCase.TestData.model); 

    % flip angles at which to test 
    thetas = [0 pi/8 pi/4 pi/2]; 
    
    % loop over flip angles computing Fisher information 
    for i=1:length(thetas)
        theta = thetas(i); 
        % use toolbox function fisher_information to compute Fisher information 
        I(i) = fisher_information(theta, testCase.TestData.model, testCase.TestData.phi);
        % use analytically-computed formula to compute Fisher information 
        I_true(i) = sin(theta)^2*testCase.TestData.phi(sin(theta)); 
    end
    
    % verify that the Fisher informations are equal 
    verifyEqual(testCase, I, I_true, 'AbsTol', epsilon())
    
    clear testCase.TestData.model
end
