function sol_dd = solve_mpc_with_dd(y_target, Sval, N, L, M_step, T, options, init_data, lambda_alpha, lambda_sigma, lambda_beta, noise)
    %#######################################################################
    % data-driven MPC for nonlinear systems, applied to CSTR
    %#######################################################################
    disp("Data-driven")
    ALPHA_PREVIOUS = true;
    sol_dd.solved = 1;
    rng(0)
    %% Setup
    % dimensions
    n = 2; % system dimsion
    nu = 2; %estimated dimension
    m = 1;
    p = 2;
    % system parameters
    h = 0.2;
    theta = 20;
    k = 300;
    M = 5;
    x_f = 0.3947;
    x_c = 0.3816;
    alpha = 0.117;
    % data parameters
    %N = 120;
    % MPC parameters
    L = L+1; % our actual horizon is L+1 shorter than this
    %M_step = 3; % number of consecutive applications of opt. input (multi-step)
    %T = 7500; % "closed-loop horizon" (simulation length)
    %x0 = [0.4;0.6];
    
    %y_T = [0.6519;0];
    y_T = [y_target; 0];
    % y_T = 0.6451;
    Q = 1;
    R = blkdiag(1,5e-2);
    S = blkdiag(Sval,0);    
    update_bound = 1e-5;
    
    % Cost matrices for problem determining alpha_sr
    %lambda_alpha_alpha = lambda_alpha;
    %lambda_sigma_alpha = lambda_sigma;
    % PE cost
    
    %% Collect initial data
    Du_cl = zeros(m,T);
    u_cl = 0.5*ones(m,T);
    u_act_cl = zeros(m,T);
    x_cl = zeros(n,T);
    y_cl = zeros(p,T);
    xi_cl = zeros(3*(m+p),T);
    
    x_cl(:,1:N+1) = init_data.x_cl(:, 1:N+1);
    u_cl(:,1:N)   = init_data.u_cl(:, 1:N);
    
    for i = 1:N
        if i>1
           Du_cl(:,i-1) = u_cl(:,i)-u_cl(:,i-1); 
        end
        y_cl(:,i) = [0 1;0 0]*x_cl(:,i)+[0;u_cl(:,i)];
        if i>=3
            u_xi = u_cl(:,i-3+1:i);
            y_xi = y_cl(:,i-3+1:i);
            xi_cl(:,i+1) = [u_xi(:);y_xi(:)]; 
        end
    end
    
    Du_cl(:,N) = -0.2+0.4*rand;
    
    u_act_cl = u_cl;
    u_cl = Du_cl;
    
    %% Set up MPC
    % Variables: alpha, u, y, us, ys, sigma, beta
    % Constraints
    u_min = 0.1;
    u_max = 2;
    us_min = 0.11;
    us_max = 1.99;
    
    
    
    % note: notation with Q and R is interchanged in the following
    % --> swap them
    Q_temp = Q;
    Q = R;
    R = Q_temp;
    
    % Cost
    H = 2*[lambda_alpha*eye(N-L+1) zeros(N-L+1,(m+p)*L+m+p+p*L+1);
           zeros(m*L,N-L+1) kron(eye(L),R) zeros(m*L,p*L) -kron(ones(L,1),R) zeros(m*L,p) zeros(m*L,p*L+1);
           zeros(p*L,N-L+1+m*L) kron(eye(L),Q) zeros(p*L,m) -kron(ones(L,1),Q) zeros(p*L,p*L+1);
           zeros(m,N-L+1) -kron(ones(1,L),R) zeros(m,p*L) L*R zeros(m,p+p*L+1);
           zeros(p,N-L+1+m*L) -kron(ones(1,L),Q) zeros(p,m) L*Q+S zeros(p,p*L+1);
           zeros(p*L,N-L+1+(m+p)*L+m+p) lambda_sigma*eye(p*L) zeros(p*L,1);
           zeros(1,N-L+1+(m+p)*L+m+p+p*L+1)];
    % penalize alpha w.r.t. zero
    f = [zeros(N-L+1+(m+p)*L+m,1);-2*S*y_T;zeros(p*L,1);-lambda_beta];
    
    % Input Constraints 
    lower = [-inf*ones(N-L+1,1);-inf*ones(m*L,1);-inf*ones(p*nu,1);repmat([-inf;u_min],L-nu,1);-inf*ones(m,1);[-inf;us_min];repmat([-inf;0],L,1);-inf];
    upper = [inf*ones(N-L+1,1);inf*ones(m*L,1);inf*ones(p*nu,1);repmat([inf;u_max],L-nu,1);inf*ones(m,1);[inf;us_max];repmat([inf;0],L,1);inf];
    
    A_ineq = [];
    b_ineq = [];
    
    fval_store = zeros(1,T);
    sol_store = zeros(N-L+1+(m+p)*L+m+p+p*L+1,T);
    time_store = zeros(1,T);
    alpha_ol = zeros(N-L+1,T);
    u_ol = zeros(m*L,T);
    y_ol = zeros(p*L,T);
    us_ol = zeros(m,T);
    ys_ol = zeros(p,T);
    sigma_ol = zeros(p*L,T);
    beta_ol = zeros(1,T);
    alpha_sr_ol = zeros(N-L+1,T);
    
    UPDATE_store = zeros(1,T);
    
    c_pe_store = zeros(1,T);
    H_ux_store = zeros(m*L+3+1,N-L+1,T);
    
    theta_store = zeros(1,T);
    
    t_update_store = [];
    % MPC iterations
    myUPDATE = true;
    for j = N+1:M_step:T-n
    

    
    
        % Cost
        % penalize alpha w.r.t. zero
        f = [zeros(N-L+1+(m+p)*L+m,1);-2*S*y_T;zeros(p*L,1);-lambda_beta];
        
        % Input Constraints 
        lower = [-inf*ones(N-L+1,1);-inf*ones(m*L,1);-inf*ones(p*nu,1);repmat([-inf;u_min],L-nu,1);-inf*ones(m,1);[-inf;us_min];repmat([-inf;0],L,1);-inf];
        upper = [inf*ones(N-L+1,1);inf*ones(m*L,1);inf*ones(p*nu,1);repmat([inf;u_max],L-nu,1);inf*ones(m,1);[inf;us_max];repmat([inf;0],L,1);inf];


    
    
    
        %j
        %% Update Constraints
        % Update data for prediction
        if j>N+1
    %       UPDATE = norm(x_cl(2,j)-y_T)>update_bound;
            UPDATE = (u_val-kron(ones(L,1),us_val))'*kron(eye(L),R)*(u_val-kron(ones(L,1),us_val))+(y_val-kron(ones(L,1),ys_val))'*kron(eye(L),Q)*(y_val-kron(ones(L,1),ys_val))>update_bound;
        else
            UPDATE = 1;
        end
        
        if UPDATE
    %     if norm(x_cl(2,j)-y_T)>0.02*norm(y_T)
            t_update_store = [t_update_store j];
            u_data = u_cl(:,j-N:j-1);
            y_data = y_cl(:,j-N:j-1);
            Hu = hankel_r(u_data(:),L,N-L+1,m);
            Hy = hankel_r(y_data(:),L,N-L+1,p);
        else
    %         UPDATE_SAVE = 0;
        end
            
        %#################################################################
            u_lin = u_act_cl(j-M_step);
            for i = 1:M_step
                u_lin = u_lin + u_cl(:,j-i);
            end
        x_curr = [x_cl(:,j);u_lin];

        if myUPDATE
            % Dynamics
            A_dyn = [Hu -eye(m*L) zeros(m*L,p*L+m+p+p*L+1);
                     Hy zeros(p*L,m*L) -eye(p*L) zeros(p*L,m+p) -eye(p*L) zeros(p*L,1)];
            b_dyn = zeros((m+p)*L,1);
            % Initial conditions
            A_init = [zeros(nu*m,N-L+1) eye(nu*m) zeros(nu*m,m*(L-nu)+p*L+m+p+p*L+1);
                      zeros(nu*p,N-L+1) zeros(nu*p,m*L) eye(nu*p) zeros(nu*p,p*(L-nu)+m+p+p*L+1)];
            u_init = u_cl(:,j-nu:j-1);
            y_init = y_cl(:,j-nu:j-1);
            b_init = [u_init(:);y_init(:)];
            % TEC
            A_TEC = [zeros((nu+1)*m,N-L+1) zeros((nu+1)*m,m*(L-nu-1)) eye((nu+1)*m) zeros((nu+1)*m,p*L) -repmat(eye(m),nu+1,1) zeros(m*(nu+1),p+p*L+1);
                     zeros((nu+1)*p,N-L+1) zeros((nu+1)*p,m*L) zeros((nu+1)*p,p*(L-nu-1)) eye((nu+1)*p) zeros((nu+1)*p,m) -repmat(eye(p),nu+1,1) zeros(p*(nu+1),p*L+1)];
            b_TEC = zeros((m+p)*(nu+1),1);
            % alpha sums to 1
            A_alpha1 = [ones(1,N-L+1) zeros(1,(m+p)*L+m+p+p*L+1)];
            b_alpha1 = 1;
            % Combine
            A_eq = [A_dyn;A_init;A_TEC;A_alpha1];
            b_eq = [b_dyn;b_init;b_TEC;b_alpha1];
            
            
            % Take alpha from last time step
            if ALPHA_PREVIOUS && j>N+1
                alpha_sr_ol(:,j) = alpha_val;
                f = [-2*lambda_alpha*alpha_sr_ol(:,j);zeros((m+p)*L+m,1);-2*S*y_T;zeros(p*L,1);-lambda_beta];
            end
            
            %% Inequality constraint for PE
            if j>N+1
                u_cand = [u_val(M_step*m+1:end);repmat(us_val,M_step,1)];
            else
                u_cand = zeros(m*L,1);
            end
            A_ineq = [zeros(1,N-L+1) 2*u_cand' zeros(1,p*L+m+p+p*L) 1];
            b_ineq = u_cand'*u_cand;
        %     A_ineq = [];
        %     b_ineq = [];
        end
        %% Solve
        [sol,fval,exitflag,output,lambda] = quadprog(H,f,A_ineq,b_ineq,A_eq,b_eq,lower,upper,[],options);
        if exitflag <=0 || isempty(sol)
            sol_dd.solved = 0;
            warning('optimization problem not solved exactly..');
            break;
        end
        
        alpha_val = sol(1:N-L+1);
        u_val = sol(N-L+1+1:N-L+1+m*L);
        y_val = sol(N-L+1+m*L+1:N-L+1+m*L+p*L);
        us_val = sol(N-L+1+(m+p)*L+1:N-L+1+(m+p)*L+m);
        ys_val = sol(N-L+1+(m+p)*L+m+1:N-L+1+(m+p)*L+m+p);
        sigma_val = sol(N-L+1+(m+p)*L+m+p+1:N-L+1+(m+p)*L+m+p+p*L);
        beta_val = sol(end);
        
        sol_store(:,j) = sol;
        fval_store(j)=fval+y_T'*S*y_T;
        alpha_ol(:,j) = alpha_val;
        u_ol(:,j) = u_val;
        y_ol(:,j) = y_val;
        us_ol(:,j) = us_val;
        ys_ol(:,j) = ys_val;
        sigma_ol(:,j) = sigma_val;
        beta_ol(:,j) = beta_val;
            
        % Simulate closed loop
        for i = j:j+M_step-1
            u_cl(:,i) = u_val(m*(nu+i-j)+1:m*(nu+i-j+1));
            u_act_cl(:,i) = u_act_cl(:,i-1)+u_cl(:,i-1);
            x_cl(:,i+1) = [x_cl(1,i)+(h/theta)*(1-x_cl(1,i))-h*k*x_cl(1,i)*exp(-M/x_cl(2,i));
                           x_cl(2,i)+(h/theta)*(x_f-x_cl(2,i))+h*k*x_cl(1,i)*exp(-M/x_cl(2,i))-h*alpha*u_act_cl(:,i)*(x_cl(2,i)-x_c)];
            if noise.on
                x_cl(:,i+1) = x_cl(:,i+1) + noise.x_cl(:,i+1); 
            end
            y_cl(:,i) = [0 1;0 0]*x_cl(:,i)+[0;u_act_cl(:,i)];
            
            u_xi = u_cl(:,i-3+1:i);
            y_xi = y_cl(:,i-3+1:i);
            xi_cl(:,i+1) = [u_xi(:);y_xi(:)]; 
        end


        if(norm(x_cl(:,j+1:j+M_step) - x_cl(:, j:j + M_step - 1)) < 1e-6) % if no change visible 
            myUPDATE = false;
        else
            myUPDATE = true;
        end
    end

    sol_dd.x_cl = x_cl;
    sol_dd.u_cl = u_act_cl;
    sol_dd.us_ol = us_ol;
    sol_dd.ys_ol = ys_ol; 
    sol_dd.fval_store = fval;
    sol_dd.sol_full = sol;
    sol_dd.time_store = time_store;
    sol_dd.params.lambda_alpha = lambda_alpha;
    sol_dd.params.lambda_sigma = lambda_sigma;
    sol_dd.params.lambda_beta = lambda_beta;
    sol_dd.params.N = N;
    sol_dd.params.L = L;
    sol_dd.y_T = y_T;
    sol_dd.noise = noise;
    sol_dd.S = S;

end