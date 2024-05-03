%% LSQ MATRIX

%% load solutions
sols_lsq = {};
for i = 1:3
    tmp = load("./simulations/sols_lsq_" + num2str(i));
    sols_lsq = [sols_lsq, tmp.solutions];
end


%% Calculate the error values

lambda_list = [0, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7]; 
N_list = [30:10:300];

dep_matrix_error = NaN*ones(length(N_list), length(lambda_list));
dep_matrix_error_conv = NaN*ones(length(N_list), length(lambda_list));

sols = sols_lsq;
sols_conv = {};

for i = 1:length(sols)
    sol = sols{i}; 
        if sol.solved
            if(sol.convergence.converged)
                dep_matrix_error_conv(find([sol.params.N == N_list]), find([sol.params.lambda == lambda_list])) = norm(sol.error(1:end-3));
            end
            dep_matrix_error(find([sol.params.N == N_list]), find([sol.params.lambda == lambda_list])) = norm(sol.error(1:end-3));

        end
end
%% Plot
close all 
[X,Y] = meshgrid(1:length(lambda_list),N_list);

pcolor(X,Y, dep_matrix_error)
%xscale log
xticklabels(split(num2str(lambda_list), '       '))
colorbar
xlabel("Parameter $\lambda$")
ylabel("Parameter $N$")
title("Error $||y - y_T||_2$: Identification-Based Algorithm")

matlab2tikz("lsqMatrix.tex")



