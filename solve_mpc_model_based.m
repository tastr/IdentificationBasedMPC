function sol_model =...
    solve_mpc_model_based(x0, target_states, H, f_cost_func, sys, S, T, M_step, L, lowerConstraint, upperConstraint, options, noise)
    disp("Model-based")
    y_T = target_states(1);
    Alin = sys.A;
    B = sys.B;
    elin = sys.e;
    Clin = sys.C;
    D = sys.D; 
    rlin = sys.r;
    f0 = sys.f0;
    g0 = sys.g0;
    sol_model.solved = true;
    
    %% initialize some variables
    [m,p] = size(D); 
    [n,~]= size(B);
    u_cl = zeros(m,T);
    x_cl = zeros(n,T);
    % initial state
    x_cl(:,1) = x0;
    
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
    
    yidx = 1; 
    for j = 1:M_step:T-n
        %j
        x_curr = x_cl(:,j);
        r_curr = zeros(p,1); 
        t1 = tic;
    
        % Linearization-based system dynamics
        A = Alin(x_curr);
        e = elin(x_curr);
        C = Clin(x_curr);
        r = rlin(x_curr);

% 
%         A = [1-h*k*exp(-M/x_curr(2))-h/theta -(M*h*k*x_curr(1)*exp(-M/x_curr(2)))/x_curr(2)^2 0;
%          h*k*exp(-M/x_curr(2)) M*h*k*x_curr(1)*exp(-M/x_curr(2))/x_curr(2)^2-h/theta+1-h*alpha*x_curr(3) -h*alpha*(x_curr(2)-x_c);
%          0 0 1];
%         B = [0;0;1];
%         e = [x_curr(1)+(h/theta)*(1-x_curr(1))-h*k*x_curr(1)*exp(-M/x_curr(2));
%             x_curr(2)+(h/theta)*(x_f-x_curr(2))+h*k*x_curr(1)*exp(-M/x_curr(2))-h*alpha*x_curr(3)*(x_curr(2)-x_c);x_curr(3)]-A*x_curr;
    
            
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
       
        %%
        time_store_jac(j) = toc(t1);
        
        t2 = tic;
        %% Solve
        [sol,fval,exitflag,output,lambda] =...
            quadprog(H,f_cost_func,A_ineq,b_ineq,A_eq,b_eq,lowerConstraint,upperConstraint,[],options);
        %fval 
        time_store_opt(j) = toc(t2);
        if exitflag <=0
            sol_model.solved = false;
            warning('optimization problem not solved exactly..') 
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
                x_cl(:,i+1) = f0(x_cl(:,i)) + B*u_cl(:,i) + [normrnd(0,noise.sigma1); normrnd(0,noise.sigma2);0]; 
            else
                x_cl(:,i+1) = f0(x_cl(:,i)) + B*u_cl(:,i); 
            end
            
        end
       
    end
    sol_model.x_cl = x_cl; 
    sol_model.u_cl = u_cl; 
    sol_model.fval = fval_store + S*y_T^2; 
    sol_model.sol = sol_store;
    sol_model.time = time_store;
    sol_model.xs_ol = xs_ol;
    sol_model.us_ol = us_ol;
    sol_model.ys_ol = ys_ol;
    sol_model.params.L = L;
    sol_model.S = S;
end