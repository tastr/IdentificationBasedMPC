%% Load data
sols_sigma = {};
for i = 1:3
    tmp = load("./simulations/sols_sigma_dd_" + num2str(i));
    sols_sigma = [sols_sigma, tmp.solutions]; 
end

sols_alpha = {};
for i = 1:3
    tmp = load("./simulations/sols_alpha_dd_" + num2str(i));
    sols_alpha = [sols_alpha, tmp.solutions]; 
end

sols_N = {};
for i = 1:2
    tmp = load("./simulations/sols_N_dd_" + num2str(i));
    sols_N = [sols_N, tmp.solutions]; 
end

%% Error for fixed lambda_sigma

lambda_alpha_list = [.5e-6:.5e-6:5e-6];
lambda_sigma_list =2.5e7;
N_list = [45, 50:10:350];%:10:200];%300:10:350;%[30, 50, 70, 90];
[dep_matrix_error_sigma, dep_matrix_error_sigma_conv, ~, min_sol_sigma] = calc_error_matrix(sols_sigma, lambda_sigma_list, lambda_alpha_list, N_list);
dep_matrix_error_sigma(:,end+1) = zeros(1,length(N_list)); % add an empty column bc matlab does not plot it otherways....
dep_matrix_error_sigma_conv(:,end+1) = zeros(1,length(N_list));

% For plots
[X_sigma,Y_sigma] = meshgrid(1:length(lambda_alpha_list)+1,N_list);
dep_matrix_error_sigma = reshape(dep_matrix_error_sigma, size(X_sigma));

xticks_sigma= 1:length(lambda_alpha_list)+1;
xticklabels_sigma=split(num2str(lambda_alpha_list))';
xticklabels_sigma{end+1}='';

%% Error for fixed lambda_alpha

lambda_alpha_list =  3.5e-6;
lambda_sigma_list = [0.01e7, 0.05e7 .1e7, .5e7, 1e7, 1.5e7, 2e7,2.5e7, 3e7 3.5e7];
N_list = [45, 50:10:350];%:10:200];%300:10:350;%[30, 50, 70, 90];

[dep_matrix_error_alpha, dep_matrix_error_alpha_conv, ~, min_sol_alpha] = calc_error_matrix(sols_alpha, lambda_sigma_list, lambda_alpha_list, N_list);

% For plots
[X_alpha,Y_alpha] = meshgrid(1:length(lambda_sigma_list)+1,N_list);
dep_matrix_error_alpha(:,1,end+1) = zeros(1,length(N_list));
dep_matrix_error_alpha_conv(:,1,end+1) = zeros(1,length(N_list));
dep_matrix_error_alpha = reshape(dep_matrix_error_alpha, size(X_alpha));

xticks_alpha= 1:length(lambda_sigma_list)+1;
xticklabels_alpha=split(num2str(lambda_sigma_list/1e7))';
for i = 1:length(xticklabels_alpha)
    xticklabels_alpha{i} = xticklabels_alpha{i} + "$\cdot 10^7$"
end
xticklabels_alpha{end+1} = '';
%% Error for fixed N

lambda_alpha_list = [.5e-6:.5e-6:5e-6];
lambda_sigma_list =[0.01e7, 0.05e7 .1e7, .5e7, 1e7, 1.5e7, 2e7,2.5e7, 3e7 3.5e7];
N_list = 120;
[X_N,Y_N] = meshgrid(1:length(lambda_sigma_list)+1,1:length(lambda_alpha_list)+1);
[dep_matrix_error_N, dep_matrix_error_N_conv, ~, min_sol_N] = calc_error_matrix(sols_N, lambda_sigma_list, lambda_alpha_list, sols_N{1}.params.N)


dep_matrix_error_N(1,:,end+1) = zeros(1,length(lambda_sigma_list));
dep_matrix_error_N(1,end+1,:) = zeros(1,length(lambda_alpha_list)+1);

dep_matrix_error_N = reshape(dep_matrix_error_N, size(X_N));

xticks_N= 1:length(lambda_sigma_list) +1;
xticklabels_N=split(num2str(lambda_sigma_list/1e7))';
for i = 1:length(xticklabels_N)
    xticklabels_N{i} = xticklabels_N{i} + "$\cdot 10^7$"
end
xticklabels_N{end+1} = '';


yticks_N= 1:length(lambda_alpha_list)+1;
yticklabels_N=split(num2str(lambda_alpha_list))';
yticklabels_N{end+1} = ''; 
%% Plot
cLowerLim = min([norm(min_sol_sigma.error), norm(min_sol_alpha.error), norm(min_sol_N.error)])
subplot(2,1,1)
pcolor(X_sigma,Y_sigma,dep_matrix_error_sigma)
colorbar
clim([cLowerLim,20])

xlabel("Parameter $\lambda_\alpha$", 'fontsize',20)
ylabel("Parameter $N$", 'fontsize',20)
title("$\lambda_\sigma = $" + num2str(sols_sigma{1}.params.lambda_sigma/1e7)+ "$\cdot 10^7$", 'fontsize',20)
xticks(xticks_sigma)
xticklabels(xticklabels_sigma)
ylim([45,300])


subplot(2,1,2)
pcolor(X_alpha,Y_alpha,dep_matrix_error_alpha)
colorbar
clim([cLowerLim,20])


xlabel("Parameter $\lambda_\alpha$", 'fontsize',20)
ylabel("Parameter $N$", 'fontsize',20)
title("$\lambda_\alpha = $" + num2str(sols_alpha{1}.params.lambda_alpha), 'fontsize',20)
xticks(xticks_alpha)
xticklabels(xticklabels_alpha)
ylim([45,300])

% 
% subplot(3,1,3)
% pcolor(X_N,Y_N,dep_matrix_error_N)
% colorbar
% clim([cLowerLim,20])
% xticks(yticks_N)
% xticklabels(yticklabels_N)
% yticks(xticks_N)
% yticklabels(xticklabels_N)
% 
% xlabel("Parameter $\lambda_\sigma$", 'fontsize',20)
% ylabel("Parameter  $\lambda_\alpha$", 'fontsize',20)
% title("$N = $" + num2str(sols_N{1}.params.N), 'fontsize',20)
% %ylim([45,300])


matlab2tikz('ddMatrix.tex')
%% Local functions 

function [dep_matrix_error, dep_matrix_error_conv, sols_conv, min_sol] = calc_error_matrix(solutions, lambda_sigma_list, lambda_alpha_list, N_list)

    dep_matrix_error =  NaN*ones(length(N_list), length(lambda_alpha_list), length(lambda_sigma_list));
    dep_matrix_error_conv =  NaN*ones(length(N_list), length(lambda_alpha_list), length(lambda_sigma_list));


    sols=solutions; 
    sols_conv = {};
    min_sol = sols{1};
    
    for i = 1:length(sols)
        sol = sols{i}; 
            if sol.solved
                if(sol.convergence.converged)
                    dep_matrix_error_conv(find([sol.params.N == N_list]), find([sol.params.lambda_alpha == lambda_alpha_list]), ...
                    find([sol.params.lambda_sigma == lambda_sigma_list])) = norm(sol.error(1:end-3));
                    sols_conv{end+1} = sol;
                end
    
                dep_matrix_error(find([sol.params.N == N_list]), find([sol.params.lambda_alpha == lambda_alpha_list]), ...
                    find([sol.params.lambda_sigma == lambda_sigma_list])) = norm(sol.error(1:end-3));
    
            end
    
    
            if (norm(sol.error)<norm(min_sol.error))
                min_sol = sol;
            end
    end

end


