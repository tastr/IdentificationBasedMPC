%#######################################################################
% Problem from BerKoeh2021_DataDrivenModelPredictiveControl
% Model-based, identification-ased and data-driven model predictive 
% control: closed-loop
% guarantees and experimental results
%#######################################################################
%% Set options & general parameters
options = mskoptimset('OptimalityTolerance', 1e-6,...
'MaxIterations', 50000,...
'ConstraintTolerance', 1e-6);%,'Display','iter');
rng('default') 
convergenceTolerance = .1; % when the solutions is assumed to converge
convergenceLength = 500; % how long the solution must be hold at y_T to be 
                         % seen as converged 

%% Where to save data
fpath = "./simulations/";
fname_base = "sols_alpha_";
%% Define Hyperparameters 
lambda_list =[];%[0,1e-14,1e-13, 1e-12];
lambda_alpha_list =  3.5e-6;%[.5e-6:.5e-6:5e-6];
lambda_sigma_list =[0.01e7, 0.05e7 .1e7, .5e7, 1e7, 1.5e7, 2e7,2.5e7, 3e7 3.5e7];
L_list = [41];
N_list = [45, 50:10:300];
solNb_mdl=1; % to save solutions
solNb_lsq=1; % to save solutions
solNb_dd=1; % to save solutions
T = 5000; % "closed-loop horizon" (simulation length)
s_bas = 100; 
y_T = 0.6519;
%% Define Noise
noise.on = false;
noise.sigma1 = 1e-6;
noise.sigma2 = 1e-6;


%% run simulation script 
simulation_main; 
