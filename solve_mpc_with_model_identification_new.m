function sol_lsq = solve_mpc_with_model_identification(x0, yrs, H,...
    f_cost_func, sys, S, T, M_step,L, lower, upper, options, Du_ident, N, init_data, reg_lsq, noise)
    disp("Identification-based")
    y_T = yrs(1); 
    A_lin = sys.A;
    B_tr = sys.B;
    e_lin = sys.e;
    C_lin = sys.C;
    D_tr = sys.D; 
    r_lin = sys.r;
    f0 = sys.f0;
    g0 = sys.g0; 

    sol_lsq.solved = 1;
    %% initialize some variables
    [m,p] = size(D_tr); 
    [n,~]= size(B_tr);
    u_cl = zeros(m,T);
    x_cl = zeros(n,T);
    y_cl = zeros(p,T);
    % initial state
    x_cl(:,1:N+1) = init_data.x_cl(:,1:N+1);
    u_cl(:,1:N) = init_data.u_cl(:,1:N);
    y_cl(:,1:N) = init_data.y_cl(:,1:N);
    
    fval_store = zeros(1,T);
    sol_store = zeros((n+m)*L+n+m+p,T);
    time_store = zeros(1,T);
    x_ol = zeros(n*L,T);
    u_ol = zeros(m*L,T);
    xs_ol = zeros(n,T);
    us_ol = zeros(m,T);
    ys_ol = zeros(p,T);
    
    A_ineq = [];
    b_ineq = [];
    UPDATE = true;

    %% Identify the matrices 
    %% Model identification 
    condNumb_state = zeros(1,T); 
    condNumb_output = zeros(1,T); 
    
    Adiff = zeros(1,T); 
    Bdiff = zeros(1,T); 

    n_last_error = zeros(1, T);

%% Solve MPC 
yidx = 1; 
    for j = N+1:M_step:T-n
        %j
        
        x_curr = x_cl(:,j);
        r_curr = zeros(p,1); 
        t1 = tic;

        if j > N && UPDATE
            [A, B, e, C, D, r, condNumb_state(j), condNumb_output(j), diff]= identify_model(x_cl(:,j-N:j), u_cl(:,j-N:j-1), y_cl(:,j-N:j-1));
            %disp("n-last-error:")
            err = 0;
            for i = j-n:j
                err = err + (norm(x_curr - x_cl(:,i)));
            end
            n_last_error(j) = err;
            Adiff(j) = diff.Adiff; 
            Bdiff(j) = diff.Bdiff; 
            Cdiff(j) = diff.Cdiff; 
            Ddiff(j) = diff.Ddiff; 
            ediff(j) = diff.ediff; 
            rdiff(j) = diff.rdiff; 
            %diff

            %% Update constraints
            % dynamics
            A_mat = [kron(eye(L-1),A) zeros((L-1)*n,n)]+[zeros((L-1)*n,n) -eye((L-1)*n)];
            A_eq_dyn = [A_mat kron(eye(L-1),B) zeros(n*(L-1),m+n+m+p)];
            b_eq_dyn = repmat(-e,L-1,1); 
            % steady-state
                % (x_s,u_s) \in Z_s 
            A_eq_ss1 = [zeros(n,(n+m)*L) eye(n)-A -B zeros(n,p)]; 
            b_eq_ss1 = e;
                % y_s = CxS + DuS + r
            A_eq_ss2 = [zeros(p,(n+m)*L) C D -eye(p)];
            b_eq_ss2 = zeros(p,1);
            % initial condition
            A_eq_init = [eye(n) zeros(n,n*(L-1)+m*L+n+m+p)];
            b_eq_init = x_curr;
            % TEC
                % - xN = x_S 
            A_eq_tec = [zeros(n,n*(L-1)) eye(n) zeros(n,m*L) -eye(n) zeros(n,m+p)];
            b_eq_tec = zeros(n,1); 
            % Combine
            A_eq = [A_eq_dyn;A_eq_ss1;A_eq_ss2;A_eq_init;A_eq_tec];
            b_eq = [b_eq_dyn;b_eq_ss1;b_eq_ss2;b_eq_init;b_eq_tec];
        end
      
        %%
                %% Solve
        [sol,fval,exitflag,output,lambda] =...
            quadprog(H,f_cost_func,A_ineq,b_ineq,A_eq,b_eq,lower,upper,[],options);
        %fval
        if exitflag <=0 || isempty(sol)
            sol_lsq.solved = 0;
            warning('optimization problem not solved exactly..') 
            if (exitflag <= 0)
                disp("Exitflag")
            end
            break
        end
            
        x_val = sol(1:n*L);
        u = sol(n*L+1:(n+m)*L);
        xs_val = sol((n+m)*L+1:(n+m)*L+n);
        us_val = sol((n+m)*L+n+1:(n+m)*L+n+m);
        ys_val = sol((n+m)*L+n+m+1:(n+m)*L+n+m+p);
        
        sol_store(:,j) = sol;
        fval_store(j)=fval+y_T'*S*y_T;
        x_ol(:,j) = x_val;
        u_ol(:,j) = u;
        xs_ol(:,j) = xs_val;
        us_ol(:,j) = us_val;
        ys_ol(:,j) = ys_val;
        
            
        % Simulate closed loop
        for i = j:j+M_step-1
            u_cl(:,i) = u(i-j+1:i-j+m);
            if noise.on
                x_cl(:,i+1) = f0(x_cl(:,i)) + B*u_cl(:,i) + [noise.x_cl(:,i); 0];
                y_cl(:,i) = g0(x_cl(:,i)) + D*u_cl(:,i) + noise.y_cl(i);
            else
                x_cl(:,i+1) = f0(x_cl(:,i)) + B*u_cl(:,i); 
                y_cl(:,i) = g0(x_cl(:,i)) + D*u_cl(:,i);
            end
        end

        % if(norm(x_cl(2,j+M_step)-y_T)<0.01*norm(y_T))
        if(norm(x_cl(:,j+1:j+M_step) - x_cl(:, j:j + M_step - 1)) < 5e-6) % if no change visible 
            UPDATE = false;
            sol_lsq.solved = j+M_step;
        else
            UPDATE = true;
        end

    end
    
    sol_lsq.x_cl = x_cl;
    sol_lsq.u_cl = u_cl;
    sol_lsq.xs_ol = xs_ol;
    sol_lsq.us_ol = us_ol;
    sol_lsq.ys_ol = ys_ol;
    sol_lsq.fval_store = fval_store;
    sol_lsq.sol_store = sol_store;
    sol_lsq.time_store  = time_store;
    sol_lsq.n_last_error = n_last_error;
    sol_lsq.condNumb_state = condNumb_state;
    sol_lsq.condNumb_output = condNumb_output; 
    sol_lsq.diff.Adiff = Adiff;
    sol_lsq.diff.Bdiff = Bdiff;
    sol_lsq.diff.Cdiff = Cdiff;
    sol_lsq.diff.Ddiff = Ddiff;
    sol_lsq.diff.ediff = ediff;
    sol_lsq.diff.rdiff = rdiff;
    sol_lsq.params.N = N;
    sol_lsq.params.L = L;
    sol_lsq.params.lambda = reg_lsq;
    sol_lsq.params.S = S;
    sol_lsq.noise = noise;
    sol_lsq.S = S;



    %% Local functions 
    function [A_lsq, B_lsq, e_lsq, C_lsq, D_lsq, r_lsq, condNumb_state, condNumb_output]= identify_model_first_time()
        
        x_ident = zeros(n,Nident+1);
        y_ident = zeros(n, Nident);        
        x_ident(:,1) =  x0;
            
        A_tr = A_lin(x0); 
        C_tr = C_lin(x0); 
        e_tr = e_lin(x0);
        r_tr = r_lin(x0); 
    
    
        for i = 1:Nident
            x_ident(:,i+1) = f0(x_ident(:,i)) + B_tr*Du_ident(:,i); 
            y_ident(:, i)  = g0(x_ident(:,i)) + D_tr*Du_ident(:,i) ;
        end
        
        [A_lsq,B_lsq,e_lsq,condNumb_state] = lsq_ext_reg(x_ident(:, 1:end-1), x_ident(:, 2:end), Du_ident,n,m,n,reg_lsq);
        [C_lsq,D_lsq,r_lsq,condNumb_output(1)] = lsq_ext_reg(x_ident(:,1:end-1), y_ident, Du_ident,n,m,p,reg_lsq);
           
        Adiff(1) = norm(A_lsq - A_tr);
        Bdiff(1) = norm(B_lsq - B_tr);
        Cdiff(1) = norm(C_lsq - C_tr);
        Ddiff(1) = norm(D_lsq - D_tr);
        ediff(1) = norm(e_lsq - e_tr);
        rdiff(1) = norm(r_lsq - r_tr); 
    
        
        %disp("||A - Alin|| = " + num2str(norm((A_lsq - A_tr))))
        %disp("||B - Blin|| = " + num2str(norm(C_lsq - C_tr)))
    end




    %% Local functions 
    function [A_lsq, B_lsq, e_lsq, C_lsq, D_lsq, r_lsq, condNumb_state, condNumb_output, diff]= identify_model(x, u, y)
        A_tr = A_lin(x(:, end)); 
        C_tr = C_lin(x(:, end)); 
        e_tr = e_lin(x(:, end));
        r_tr = r_lin(x(:, end)); 
    
    
    
        [A_lsq,B_lsq,e_lsq,condNumb_state] = lsq_ext_reg(x(:, 1:end-1), x(:, 2:end), u,n,m,n,reg_lsq);
        [C_lsq,D_lsq,r_lsq,condNumb_output(1)] = lsq_ext_reg(x(:,1:end-1), y, u,n,m,p,reg_lsq);
           
        diff.Adiff = norm(A_lsq - A_tr);
        diff.Bdiff = norm(B_lsq - B_tr);
        diff.Cdiff = norm(C_lsq - C_tr);
        diff.Ddiff = norm(D_lsq - D_tr);
        diff.ediff = norm(e_lsq - e_tr);
        diff.rdiff = norm(r_lsq - r_tr); 
    
        
%         disp("||A - Alin|| = " + num2str(norm((A_lsq - A_tr))))
%         disp("||B - Blin|| = " + num2str(norm(B_lsq - B_tr)))
%         disp("||C - Clin|| = " + num2str(norm((C_lsq - C_tr))))
%         disp("||D - Dlin|| = " + num2str(norm(D_lsq - D_tr)))
%         disp("||e - elin|| = " + num2str(norm((e_lsq - e_tr))))
%         disp("||r - rlin|| = " + num2str(norm(r_lsq - r_tr)))
    end
end