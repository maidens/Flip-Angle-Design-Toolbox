%SCRIPT DEMO provides a demonstration of the capabilities of the Flip Angle Design Toolbox 
%   Uses two- or three-compartment model of a simple pathway 
%
%   Flip Angle Design Toolbox 
%   John Maidens (maidens@eecs.berkeley.edu) 
%   June 2014 
 
clear all
close all
clc

% verify that required toolboxes are installed 
check_system_requirements(); 


%% Specify system model 

% initialize model object 
model = linear_exchange_model; 

% define model parameters
syms R1P R1L kPL kTRANS 
% define input parameters 
syms t0 alpha_1 beta_1 A0 
% define initial state parameters
syms P0 L0 

% parameters of interest 
% (those for which we wish to compute an estimate with minimal variance) 
model.parameters_of_interest                = [kPL  kTRANS]; 
model.parameters_of_interest_nominal_values = [0.02 0.04]; 

% nuisance parameters
% (those parameters that are unknown but whose estimates we only care about
% insofar as they allow us to estamate the parameters of interest) 
model.nuisance_parameters                = [alpha_1 beta_1 A0];
model.nuisance_parameters_nominal_values = [2       5      1 ]; 

% known parameters
% (those whose values are assumed to be known constants) 
model.known_parameters       = [R1P  R1L  t0 P0 L0 ]; 
model.known_parameter_values = [1/30 1/30  0  0  0 ];  

% define system matrices for differential eq. 
%   dx/dt = A*x(t) + B*u(t)
%    y(t) = C*x(t) + D*u(t) 

% two-site exchange model with input feedthrough 
model.A = [ -kTRANS-R1P     0      0  ;
             kTRANS     -kPL-R1P   0  ;
                0          kPL   -R1L];  
         
model.B = [1; 0; 0]; 

model.C = [1 0 0; 
           0 1 0; 
           0 0 1]; 
       
model.D = [0; 
           0;
           0]; 

% define input function shape  
model.u = @(t) A0 * (t - t0)^alpha_1 *exp(-(t - t0)/beta_1); % gamma-variate input  
% model.u = @(t) 10*rectangularPulse(0, 15, t);              % boxcar input 

% define initial condition 
% model.x0 = [P0; L0]; 
model.x0 = [0; P0; L0]; 

% define repetition time
model.TR = 2; 

% define number of acquisitions 
model.N = 25; 

% choose noise type
model.noise_type = 'Rician';
% model.noise_type = 'None';

% choose noise magnitude  
model.noise_parameters = [0.1 0.1 0.1]; % sigma^2 values for the noise 

% choose flip angle input matrix 
%   This allows you to set linear equality constraints on the flip angles
%   for example setting: 
%
%      model.flip_angle_input_matrix = [1 0; 
%                                       0 1; 
%                                       1 0]; 
%
%   fixes the first and third flip angles to be equal one another. 
%   Consider defaulting to
% 
%      model.flip_angle_input_matrix = eye(model.n) 
% 
%   if you wish to compute all flip angles separately. 
model.flip_angle_input_matrix = [1 0; 
                                 1 0;
                                 0 1];  
                             
% model.flip_angle_input_matrix = eye(model.m + model.n)                              

% choose design criterion 
design_criterion = 'D-optimal'; 
% design_criterion = 'E-optimal'; 
% design_criterion = 'A-optimal'; 
% design_criterion = 'T-optimal'; 
% design_criterion = 'totalSNR'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
model = sensitivities(model);  


%% Design optimal flip angles

% specify optimization start point and options for MATLAB optimization toolbox 
initial_q_value = 5*pi/180*ones(size(model.flip_angle_input_matrix, 2), model.N);
options = optimset('MaxFunEvals', 5000, 'MaxIter', 500, 'Display', 'iter');

% specify regularization parameter
% larger values of lambda lead to smoother flip angle sequence 
% set lambda = 0 for unregularized optimization 
lambda = 0.1; 

% perform optimization 
[thetas, ~, q_opt] = optimal_flip_angle_design_regularized(model, design_criterion, ...
    initial_q_value, lambda, options); 


%% Plot optimized flip angles 
figure 
plot(q_opt'.*180./pi, 'x-') 
title('Optimized flip angle sequence') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('Pyr', 'Lac')
% axis([1 model.N 0 100])

thetas_opt = thetas(:, 1:2); 


%% Generate simulated data from model 

% generate simulated trajecories
[y, y_true] = generate_data(model, thetas); 

% plot simulated trajectories 
figure
plot(model.TR*(0:model.N-1), y', 'o-')
title('Simulated data') 
xlabel('time (s)')
ylabel('measured magnetization (au)')
legend('Pyr', 'Lac', 'AIF')


%% Estimate model parameters from simulated data 

% choose loss function for parameter fit 
goodness_of_fit_criterion = 'maximum-likelihood'; 
% goodness_of_fit_criterion = 'least-squares'

% fit parameters values to simulated data 
[parameters_of_interest_est, nuisance_parameters_est] ...
    = parameter_estimation(y, model, goodness_of_fit_criterion, thetas) 


%% Plot results of estimate 

model.parameters_of_interest_nominal_values = parameters_of_interest_est; 
model.nuisance_parameters_nominal_values = nuisance_parameters_est; 
model = discretize(model); 
[~, y_fit] = generate_data(model, thetas); 

figure
plot(model.TR*(0:model.N-1), y', 'o-', (model.TR*(0:model.N-1)), y_fit')
title('Model fit to simulated data') 
xlabel('time (s)')
ylabel('measured magnetization (au)')
legend('Pyr data', 'Lac data', 'AIF data', 'Pyr fit', 'Lac fit', 'AIF fit')


