%SCRIPT DEMO provides a demonstration of the capabilities of the Flip Angle Design Toolbox 
%   Uses two- or three-compartment model of a simple pathway 
%
%   Flip Angle Design Toolbox 
%   John Maidens (maidens@eecs.berkeley.edu) 
%   June 2014 
 
clear all
close all
clc

% check that required toolboxes are installed 
v = ver;
toolboxes = setdiff({v.Name}, 'MATLAB'); 
if ~strncmp('Optimization Toolbox', toolboxes, 20)
    error('Flip Angle Design Toolbox requires MATLAB Optimization Toolbox')
elseif ~strncmp('Symbolic Math Toolbox', toolboxes, 21)
    error('Flip Angle Design Toolbox requires MATLAB Symbolic Math Toolbox')
end

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
model.parameters_of_interest_nominal_values = [0.05 0.04]; 

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

% three-site exhange model 
% model.A = [ -kPL-kPA-R1P  0    0;
%             kPL         -R1L  0;
%             kPA          0   -R1A];   
% model.B = [kTRANS; 0; 0]; 

% two-site exchange model 
model.A = [ -kPL-R1P  0   ;
             kPL     -R1L];   
model.B = [kTRANS; 0]; 

% define input function shape  
model.u = @(t) A0 * (t - t0)^alpha_1 *exp(-(t - t0)/beta_1); % gamma-variate input  
% model.u = @(t) 10*rectangularPulse(0, 15, t);              % boxcar input 

% define initial condition 
model.x0 = [P0; L0]; 

% define repetition time
model.TR = 2; 

% define number of acquisitions 
model.N = 25; 

% choose noise type
model.noise_type = 'Rician';
% model.noise_type = 'None';
if model.n == 3
    model.noise_parameters = [0.01 0.01 0.01 0.1]; 
elseif model.n == 2
    model.noise_parameters = [0.01 0.01 0.1]; 
end

% choose design criterion 
% design_criterion = 'totalSNR'; 
design_criterion = 'D-optimal'; 
% design_criterion = 'E-optimal'; 
% design_criterion = 'T-optimal'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
if ~model.sensitivities_computed ...
        && (strcmp(design_criterion, 'D-optimal') ...
            || strcmp(design_criterion, 'E-optimal') ...
            || strcmp(design_criterion, 'T-optimal') ...
           )
    model = sensitivities(model);  
end

% design optimal flip angles for maximum likelihood estimation
initial_thetas_value = pi/2*ones(model.N, model.n + model.m);
options = optimset('MaxFunEvals', 5000, 'MaxIter', 200, 'Display', 'iter'); 
thetas = optimal_flip_angle_design(model, design_criterion, ...
    initial_thetas_value, options); 

% plot optimal flip angles 
figure 
plot(thetas.*180./pi, 'x-') 
title('Optimal flip angle scheme') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
if model.n == 2
    legend('Pyr', 'Lac', 'AIF')
else
    legend('Pyr', 'Lac', 'Ala', 'AIF') 
end
axis([1 model.N 0 100])

% generate simulated trajecories
[y, y_true] = generate_data(model, thetas); 

% choose loss function for parameter fit 
goodness_of_fit_criterion = 'maximum-likelihood'; 
% goodness_of_fit_criterion = 'least-squares'

% fit parameters values to simulated data 
[parameters_of_interest_est, nuisance_parameters_est] ...
    = parameter_estimation(y, model, goodness_of_fit_criterion, thetas) 

% plot simulated trajectories and corresponding fit 
figure
plot(model.TR*(1:model.N), y', 'o-')
title('Simulated data') 
xlabel('time (s)')
ylabel('measured magnetization (au)')
if model.n == 2
    legend('Pyr', 'Lac', 'AIF')
else
    legend('Pyr', 'Lac', 'Ala', 'AIF') 
end


