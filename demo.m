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
syms R1P R1L R1A kPL kPA kTRANS 
% define input parameters 
syms t0 alpha_1 beta_1 A0 
% define noise parameters 
syms sigma_1 sigma_2 sigma_3
% define initial state parameters
syms P0 L0 

% parameters of interest 
% (those for which we wish to compute an estimate with minimal variance) 
model.parameters_of_interest = [kPL kTRANS]; 
model.parameters_of_interest_nominal_values = [0.02 0.04]; 

% nuisance parameters
% (those parameters that are unknown but whose estimates we only care about
% insofar as they allow us to estamate the parameters of interest) 
model.nuisance_parameters = [alpha_1 beta_1 A0];
model.nuisance_parameters_nominal_values = [ 2  5  1]; 

% known parameters
% (those whose values are assumed to be known constants) 
model.known_parameters = [R1P R1L R1A t0 P0 L0]; 
model.known_parameter_values = [1/35 1/30 1/25 0 0 0];  

% define system matrices for differential eq. dx/dt = A*x(t) + B*u(t)
% the input should be the same chemical compound as the first state 

% two-site exchange model 
model.A = [ -kPL-R1P  0   ;
             kPL     -R1L];   
model.B = [kTRANS; 0]; 

% three-site exhange model 
% model.A = [ -kPL-kPA-R1P  0    0;
%             kPL         -R1L  0;
%             kPA          0   -R1A];   
% model.B = [kTRANS; 0; 0]; 

% define input function shape  
model.u = @(t) A0 * (t - t0)^alpha_1 *exp(-(t - t0)/beta_1); % gamma-variate input  
% model.u = @(t) 10*rectangularPulse(0, 15, t);              % boxcar input 

% define initial condition 
model.x0 = [P0; L0]; 
% model.x0 = [P0; L0; 0]; 

% define repetition time
model.TR = 2; 

% define number of acquisitions 
model.N = 25; 

% choose noise type
model.noise_type = 'Rician';
% model.noise_type = 'None';

model.noise_parameters = [0.01 0.01 0.1]; 
% model.noise_parameters = [0.01 0.01 0.01 0.1]; 

% choose design criterion 
design_criterion = 'D-optimal'; 
% design_criterion = 'E-optimal'; 
% design_criterion = 'A-optimal'; 
% design_criterion = 'T-optimal'; 
% design_criterion = 'totalSNR'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
if ~model.sensitivities_computed ...
        && (strcmp(design_criterion, 'D-optimal') ...
            || strcmp(design_criterion, 'E-optimal') ...
            || strcmp(design_criterion, 'A-optimal') ...
            || strcmp(design_criterion, 'T-optimal') ...
           )
    model = sensitivities(model);  
end




%% Design optimal flip angles

% specify optimization start point and options for MATLAB optimization toolbox 
initial_thetas_value = pi/2*ones(model.N, model.n);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 200, 'Display', 'iter');

% perform optimization 
thetas = optimal_flip_angle_design(model, design_criterion, ...
    initial_thetas_value, options); 

% plot optimal flip angles 
figure 
plot(thetas.*180./pi, 'x-') 
title('Optimal flip angle scheme') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('Pyr', 'Lac', 'AIF')
axis([1 model.N 0 100])

thetas_opt = thetas(:, 1:2); 
save('flip_angles.mat', 'thetas_opt')

% 
% 
% %% Generate simulated data from model 
% 
% % generate simulated trajecories
% [y, y_true] = generate_data(model, thetas); 
% 
% % plot simulated trajectories 
% figure
% plot(model.TR*(1:model.N), y', 'o-')
% title('Simulated data') 
% xlabel('time (s)')
% ylabel('measured magnetization (au)')
% legend('Pyr', 'Lac', 'AIF')
% 
% 
% 
% 
% %% Estimate model parameters from data 
% 
% % choose loss function for parameter fit 
% goodness_of_fit_criterion = 'maximum-likelihood'; 
% % goodness_of_fit_criterion = 'least-squares'
% 
% % fit parameters values to simulated data 
% [parameters_of_interest_est, nuisance_parameters_est] ...
%     = parameter_estimation(y, model, goodness_of_fit_criterion, thetas) 
% 
% 
% 

