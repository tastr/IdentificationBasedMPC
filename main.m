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
lambda_list =0;%[0,1e-14,1e-13, 1e-12];
lambda_alpha_list = 0;%[.5e-6, 3e-6, 5]
lambda_sigma_list =0;%[.3e7, .5e7:.5e7:5e7];
L_list = [41];
N_list = 50;%[45, 50:10:350];
solNb_mdl=1; % to save solutions
solNb_lsq=1; % to save solutions
solNb_dd=1; % to save solutions
T = 5000; % "closed-loop horizon" (simulation length)
s_bas = 10; 
y_T = 0.6519;
%% Define Noise
noise.on = true;
noise.sigma1 = 1e-6;
noise.sigma2 = 1e-6;

%% MPC parameters
% system dimensions
n = 3;
m = 1;
p = 1;

M_step = 3; % number of consecutive applications of opt. input (multi-step)
Q = blkdiag(1,1,5e-2);
R = 1e0*eye(m);
S = s_bas*eye(p);

%Constraint parameters
u_min = 0.1;
u_max = 2;
us_min = 0.11;
us_max = 1.99;


%% Define system 



% system parameters
h = 0.2;
theta = 20;
k = 300;
M = 5;
x_f = 0.3947;
x_c = 0.3816;
alpha = 0.117;

% Init + Target parameters  
x0 = [0.9492;0.43;.1366];
x_eq = [0.2632;0.6519;0.7583]; % used for tests 

% system 
f0 = @(x)[x(1)+(h/theta)*(1-x(1))-h*k*x(1)*exp(-M/x(2));
            x(2)+(h/theta)*(x_f-x(2))+h*k*x(1)*exp(-M/x(2))-h*alpha*x(3)*(x(2)-x_c);
            x(3)];
g0 = @(x) x(2);

% Linearized system 
A_lin = @(x)[1-h*k*exp(-M/x(2))-h/theta -(M*h*k*x(1)*exp(-M/x(2)))/x(2)^2 0;
     h*k*exp(-M/x(2)) M*h*k*x(1)*exp(-M/x(2))/x(2)^2-h/theta+1-h*alpha*x(3) -h*alpha*(x(2)-x_c);
    0 0 1];
B = [0; 0; 1]; 
C_lin = @(x)[0 1 0];
D = 1;
e_lin = @(x)[x(1)+(h/theta)*(1-x(1))-h*k*x(1)*exp(-M/x(2));
       x(2)+(h/theta)*(x_f-x(2))+h*k*x(1)*exp(-M/x(2))-h*alpha*x(3)*(x(2)-x_c);x(3)]-A_lin(x)*x;
r_lin = @(x)0; 



sys.A =A_lin; 
sys.B = B; 
sys.C = C_lin;
sys.D = 0;
sys.e = e_lin;
sys.r = r_lin; 
sys.f0 = f0; 
sys.g0 = g0;

% predefine solutions cell-arrays
sols_dd = {};
sols_lsq = {};
sols_mdl = {};

% load noise
noiseValues = load('./datasets/initNoise.mat');% the values for rand-values
                                               % to represent noise are 
                                               % saved here to ensure 
                                               % repeatability
initNoise = noiseValues.initNoise; 
noise.x_cl = [noise.sigma1*initNoise.rand1; noise.sigma1*initNoise.rand2];
noise.y_cl = [noise.sigma1*initNoise.rand3]; % is mostly not important as 
                                             % if we know the full state we
                                             % probably can recostruct output; 

% Load data for initialization
initData = load('./datasets/initData.mat');
init_data_lsq.x_cl = initData.sol_model.x_cl(:,:);
init_data_lsq.u_cl = initData.sol_model.u_cl(:,:);
init_data_lsq.y_cl = initData.sol_model.x_cl(2,:);

init_data_dd.x_cl = init_data_lsq.x_cl(1:2,:);
init_data_dd.u_cl = init_data_lsq.x_cl(3,:);
init_data_dd.y_cl = init_data_lsq.x_cl(2,:);


%% Start simulations for different hyperparameters
% how many simulations are to be done 
totalSimNumber = length(L_list)+ ...
                 length(L_list)*length(N_list)*length(lambda_list) + ...
                 length(L_list)*length(N_list)*...
                 length(lambda_sigma_list)*length(lambda_alpha_list); 
simNumber = 0; 

for l = 1:length(L_list)
    L = L_list(l);

    % Cost
    H = 2*[kron(eye(L),Q) zeros(n*L,m*L) -kron(ones(L,1),Q) zeros(n*L,m+p);
           zeros(m*L,n*L) kron(eye(L),R) zeros(m*L,n) -kron(ones(L,1),R) zeros(m*L,p);
           -kron(ones(1,L),Q) zeros(n,m*L) L*Q zeros(n,m+p);
           zeros(m,n*L) -kron(ones(1,L),R) zeros(m,n) L*R zeros(m,p);
           zeros(p,(n+m)*L+n+m) S];
    f_cost_func = [zeros((n+m)*L+n+m,1);-2*S*y_T];
    
    % Input Constraints 
    lower = [repmat([-inf;-inf;u_min],L,1);-inf*ones(m*L,1);-inf;-inf;us_min;-inf*ones(m,1);-inf*ones(p,1)];
    upper = [repmat([inf;inf;u_max],L,1);inf*ones(m*L,1);inf;inf;us_max;inf*ones(m,1);inf*ones(p,1)];
    
    A_ineq = [];
    b_ineq = [];

    
    %% Model based solution 
    sol_model = solve_mpc_model_based(x0, y_T, H, f_cost_func, sys, S, T, M_step, L, lower, upper, options, noise);
    sols_mdl{end+1} = sol_model;
    
    simNumber=simNumber+1;
    disp("Simulation " + num2str(simNumber) + " of " + num2str(totalSimNumber))
    
    % save and release memory 
    if length(sols_dd) == 100
        saveSimulations(sols_mdl, solNb_mdl, "mdl", fpath, fname_base)
        solNb_mdl = solNb_mdl + 1;
        sols_mdl = {};
    end

    for n_idx = 1:length(N_list)
        %% Least Squares Solution
        N = N_list(n_idx);
        

        for lam_idx = 1:length(lambda_list)
            lambda = lambda_list(lam_idx);
            sol_lsq = solve_mpc_with_model_identification_new(x0, y_T, H, f_cost_func, sys, S, T, M_step,L, lower, upper, options, NaN, N, init_data_lsq, lambda, noise);
            error = sol_lsq.x_cl(2,:) - y_T;
            error_list_lsq(lam_idx) = norm(error);
            conv = all(abs(sol_lsq.x_cl(2,end-convergenceLength:end-3) - y_T)...
                       <convergenceTolerance); % end-3 because T is not multiple of M_step
            sol_lsq.error = error;
            sol_lsq.convergence.tolerance = convergenceTolerance;
            sol_lsq.convergence.converged = conv;

            sols_lsq{end+1} = sol_lsq; 
            simNumber=simNumber+1;
            disp("Simulation " + num2str(simNumber) + " of " + num2str(totalSimNumber))

            disp("Identification-Based")
            disp("N = " + num2str(N))
            disp("lambda = " + num2str(lambda))
            disp("solved: " + sol_lsq.solved)
            disp("convergence: " + conv)


            % save and release memory 
            if length(sols_lsq) == 100
                saveSimulations(sols_lsq, solNb_lsq, "lsq", fpath, fname_base)
                solNb_lsq = solNb_lsq + 1;
                sols_lsq = {};
            end


        end



        
        %% Data-driven solution

        if N <= L+n % required according to theory 
            continue;
        end
        lambda_beta = 0;
        error_matrix_dd = ones(length(lambda_alpha_list), length(lambda_sigma_list));
        for lam_al_idx = 1:length(lambda_alpha_list)
            for lam_sig_idx = 1:length(lambda_sigma_list)
                lambda_alpha = lambda_alpha_list(lam_al_idx);
                lambda_sigma = lambda_sigma_list(lam_sig_idx);
                sol_dd = solve_mpc_with_dd(y_T, S, N, L, M_step, T, options, init_data_dd, lambda_alpha, lambda_sigma, lambda_beta, noise);

                error = sol_dd.x_cl(2,:) - y_T;
                sol_dd.error = error;
                error_matrix_dd(lam_al_idx, lam_sig_idx) = norm(error); 
                conv = all(abs(sol_dd.x_cl(2,end-convergenceLength:end-3) - y_T)...
                       <convergenceTolerance); % end-3 because T is not multiple of M_step
                sol_dd.error = error;
                sol_dd.convergence.tolerance = convergenceTolerance;
                sol_dd.convergence.converged = conv;


                disp("Data-Driven")
                disp("N = " + num2str(N))
                disp("lambda_alpha = " + num2str(lambda_alpha))
                disp("lambda_sigma = " + num2str(lambda_sigma))
                disp("solved: " + sol_dd.solved)
                disp("convergence: " + conv)
                
                sols_dd{end+1} = sol_dd;
                simNumber=simNumber+1;
                disp("Simulation " + num2str(simNumber) + " of " + num2str(totalSimNumber))

                % save and release memory 
                if length(sols_dd) == 100
                    saveSimulations(sols_dd, solNb_dd, "dd", fpath, fname_base)
                    solNb_dd = solNb_dd + 1;
                    sols_dd = {};
                end
            end
        end
        end

 
end


% save the last simulations
saveSimulations(sols_mdl, solNb_mdl, "mdl", fpath, fname_base)
saveSimulations(sols_lsq, solNb_lsq, "lsq", fpath, fname_base)
saveSimulations(sols_dd, solNb_dd, "dd", fpath, fname_base)


%% %%%%%%%%%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveSimulations(solutions, solNb, type, fpath, fname_base)
    fname = fpath + fname_base + type+ "_" + num2str(solNb) + ".mat"; 
    while isfile(fname)
        fname = fname + "_1";
    end
    save(fname, "solutions");
end
