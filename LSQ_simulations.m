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
fname_base = "sols_";
%% Define Hyperparameters 
lambda_list =[0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]; 
lambda_alpha_list =  [];
lambda_sigma_list =[];
L_list = [41];
N_list = [20:5:45];%[45, 50:10:300];
solNb_mdl=0; % to save solutions
solNb_lsq=3; % to save solutions
solNb_dd=0; % to save solutions
T = 5000; % "closed-loop horizon" (simulation length)
s_bas = 100; 
y_T = 0.6519;
%% Define Noise
noise.on = false;
noise.sigma1 = 1e-6;
noise.sigma2 = 1e-6;

%% run simulation script 
simulation_main; 