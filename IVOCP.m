function[ ] = IVOCP(mpcc_mode, rho_vals, nxsamples, N, outdir)
% IVOCP     Solve the IVOCP discretization with various modes.
%           Results are stored in .mat files to the directory 'outdir'.
% Inputs:
%   mpcc_mode   : An integer representing the mode.
%                   - 1: LCQP
%                   - 2: LCQP sparse
%                   - 3: IPOPT Penalty
%                   - 4: IPOPT Smoothed
%                   - 5: IPOPT Relaxed
%   rho_vals  : A double array containing the penalty homotopy.
%   nxsamples   : An integer describing the amount of initializations to
%                   be solved.
%   N           : The number of discretization nodes.
%   outdir      : The output directory.
%

%% Prepare LCQP
% Import casadi
import casadi.*

% Time horizon and step length
T = 2;
h = T/N;

% Load CasADi collocation matrices
collocation_times;

% State variables
x = MX.sym('x');

% Algebraic variables
y = MX.sym('y');
lambda0 = MX.sym('lambda0');
l = MX.sym('l');
z = [y,lambda0,l];

% Relaxation parameter
sigma = MX.sym('sigma');

% Regularization factor
regTerm = eps;

% State/Algebraic dimensions
nx = length(x);
nz = length(z);

% State and algebraic state bounds
lb_x = -inf;
ub_x = inf;
lb_z  = [0;0;0];
ub_z = [1;inf;1];

% Constraints and objective contributions
f_x = 1*y+3*l;
f_q = x^2;
f_z = l - (1 - y);
F = Function('f_x',{x, z},{f_x, f_q, f_z});

% Initialize objective expression 
J = 0;

% Initialize variables
w = {};
w0 = [];

% Initialize (box and linear) constraints
g = {};
lbw = [];
ubw = [];
lbg = [];
ubg = [];

% Values for building initial guess 
% (Initial guess of state is adapted later)
x0 = -1;
z0 = [0; 1; 1];

% "Lift" initial conditions
Xk = MX.sym('X0', nx);
w = {w{:}, Xk};

% Add regularization term to lifted variable
J = J + regTerm*(Xk'*Xk);

% Add bounds for lifted variable
lbw = [lbw; -inf];
ubw = [ubw; inf];

% Add initial guess for lifted variable
w0 = [w0; x0];

% Keep track of indices
ind_x = 1:nx;
ind_z = [];
ind_total = ind_x;

% Total number of optimization variables and complementarity pairs
NX = nx + N*(nx + nz);
NCOMP = N*2;

% LCQP requires complementarity matrices
S1 = zeros(NCOMP, NX);
S2 = zeros(NCOMP, NX);
compCounter = 1;

% Smoothed/Relaxed strategy require adjusted algebraic bounds
if (mpcc_mode == 4)
    lbComp = ones(2,1);
    ubComp = ones(2,1);
elseif (mpcc_mode == 5)
    lbComp = -inf(2,1);
    ubComp = ones(2,1);
end

%% Formulate the NLP
for k=0:N-1
    % New node
    Xkj = MX.sym(['X_' num2str(k)], nx);
    Zkj = MX.sym(['Z_' num2str(k)], nz);
    w = {w{:}, Xkj, Zkj};

    % box constraints for new node
    lbw = [lbw; lb_x; lb_z];
    ubw = [ubw; ub_x; ub_z];
    
    % initial guess for new node
    w0 = [w0; x0; z0];

    % Update (state) indices
    newStates = ind_total(end)+1:ind_total(end)+nx; 
    ind_x = [ind_x, newStates];
    ind_total  = [ind_total, newStates];

    % Update (algebraic state) indices
    newComplementarities = ind_total(end)+1:ind_total(end)+nz;        
    ind_z = [ind_z, newComplementarities];
    ind_total  = [ind_total,newComplementarities];
    
    % Build complementarity matrices (required for easy compl. evaluation)
    % Get new indices
    yInd = newComplementarities(1);
    lam0Ind = newComplementarities(2);
    lInd = newComplementarities(3);
    xInd = newStates(1);

    % y*lambda0
    S1(compCounter, yInd) = 1;
    S2(compCounter, lam0Ind) = 1;
    compCounter = compCounter + 1;

    % (1-y)*lambda1
    S1(compCounter, lInd) = 1;
    S2(compCounter, [xInd, lam0Ind]) = [1, 1];
    compCounter = compCounter + 1;
    
    % Collocation equation
    Xk_end = D(1)*Xk;
    xp = C(1,2)*Xk + C(2,2)*Xkj;
    
    % Obtain algebraic variables
    y = Zkj(1);
    lam0 = Zkj(2);
    l = Zkj(3);
            
    % Complementarity constraints (depends on mpcc_mode)
    if (mpcc_mode == 1 || mpcc_mode == 2)
        % Lambda1 is dependent on mpcc strategy:
        lam1 = Xkj + lam0;
    elseif (mpcc_mode == 3)
        % IPOPT Penalty
        % Lambda1 is dependent on mpcc strategy:
        lam1 = Xkj + lam0;
        
        % Add penalty term to objective
        J = J + sigma*(y * lam0);
        J = J + sigma*(l * lam1);
        
    else
        %IPOPT smoothing relaxing
        lam1 = Xkj/sigma + lam0;
        
        % Add smoothing/relaxation terms to constraints
        gj = [ y*lam0; ...
               (1-y)*lam1];
        g = {g{:}, gj};
        lbg = [lbg; lbComp];
        ubg = [ubg; ubComp];     
    end
    
    % Linear constraints on algbraic vars
    % 1) Positivity of lambda1:
    g = {g{:}, lam1};
    lbg = [lbg; 0];
    ubg = [ubg; inf];
    
    % 2) In smoothing the following would lead to overconstrainedness.
    if (mpcc_mode ~= 4)
        g = {g{:}, l - (1 - y)};
        lbg = [lbg; 0];
        ubg = [ubg; 0];
    end

    % ODE and objective
    [fx, fq, fz] = F(Xkj, Zkj);

    % collect discretized DAE equations
    g = {g{:}, h*fx - xp, fz};
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];

    % Add contribution to the end state
    Xk_end = Xk_end + D(2)*Xkj;

    % Add contribution to quadrature function
    J = J + B(2)*fq*h;

    % Add regularization term to algebraic variables
    J = J + regTerm*(Zkj'*Zkj);

    % Keep track of this node for next node
    Xk = Xk_end;
end

% Terminal cost
J = J + (Xk_end-5/3)^2;

% Create array containing all optimization variables
x = vertcat(w{:});

% Original objective function
J_fun = Function('J_fun', {x, sigma}, {J});
fprintf('Peparing the NLP done ... \n');

%% Initialize the solvers
% Constraint matrix
constr = vertcat(g{:});
if (mpcc_mode == 1 || mpcc_mode == 2)
    % Create QP from NLP
    Constr = Function('Constr', {x}, {constr});
    A_Fun = Function('A_Fun', {x}, {jacobian(constr, x)});
    A = full(A_Fun(zeros(size(x))));
    
    % Linearization correction term
    constr_constant = Constr(zeros(size(x)));
    lbg = lbg - full(constr_constant);
    ubg = ubg - full(constr_constant);
    
    % Linear objective term
    J_Jac_Fun = Function('J_Jac_fun', {x}, {jacobian(J, x)});
    g = full(J_Jac_Fun(zeros(size(x))))';

    % Quadratic objective term
    Q_Fun = Function('Q_fun', {x}, {hessian(J, x)});
    Q = full(Q_Fun(zeros(size(x))));
        
    % LCQP parameters
    params.x0 = w0;
    params.initialPenaltyParameter = rho_vals(2);
    params.solveZeroPenaltyFirst = true;
    params.penaltyUpdateFactor = rho_vals(3)/rho_vals(2);
    params.maxIterations = 1000;
    params.printLevel = 0;
    params.Rbreak = rho_vals(end);
    
    if (mpcc_mode == 2)
        params.useSparseMatrices = true;
    end
else
    % Create IPOPT solver through CasADi
    opts_ipopt.verbose = false;
    opts_ipopt.ipopt.print_level = 0;
    opts_ipopt.print_time = 0;
    opts_ipopt.print_out = 0;
    prob = struct('f', J, 'x', x, 'g', constr, 'p', sigma);
end

% Solution vectors 
w_opt_vec = zeros(NX, nxsamples);
compl_vals = zeros(nxsamples, 1);
obj_originals = zeros(nxsamples, 1);
time_vals = zeros(nxsamples, 1);

%% Solution sequence
% Linspace for fixed x0 values
xx = linspace(-1.9,-0.9,nxsamples)';
C = [S1; S2]'*[S2; S1];

for i = 1:nxsamples 
    % Obtain initial guess for states from different forward simulations
    initGuess = BuildInitialGuessEuler(@(x) 2 - sign(x), xx(i), h, N+1);    
    w0(ind_x) = initGuess;
    params.x0 = w0;
            
    if (mpcc_mode == 1 || mpcc_mode == 2)
        % Solve using LCQP
        tic;
        w_opt = LCQP(Q, g, S1, S2, A, lbg, ubg, lbw, ubw, params);
        time_vals(i) = toc;
    else
        % Solve using IPOPT
        
        % Reset solver
        prob = struct('f', J, 'x', x, 'g', constr, 'p', sigma);
        solver = nlpsol('solver', 'ipopt', prob, opts_ipopt);
        
        % Go through homotopy
        for j = 1:length(rho_vals)        
            sigma_val = rho_vals(j);
            
            tic;
            sol = solver('x0', params.x0, 'lbx', lbw, 'ubx', ubw,...
                         'lbg', lbg, 'ubg', ubg, 'p', sigma_val);
            time_vals(i) = time_vals(i) + toc;  
            w_opt = full(sol.x);  
            
            params.x0 = w_opt;
        end
    end
    
    % Save variables
    w_opt = full(w_opt);
    obj_originals(i) = full(J_fun(w_opt, 0));
    w_opt_vec(:, i) = w_opt;
    compl_vals(i) = w_opt'*S1'*S2*w_opt;
end

%% Save results to mat file
if (mpcc_mode == 1)
    algo = "LCQP";
elseif (mpcc_mode == 2)
    algo = "LCQP_Sparse";
elseif (mpcc_mode == 3)
    algo = "nlp_pen";
elseif (mpcc_mode == 4)
    algo = "smoothed";
elseif (mpcc_mode == 5)
    algo = "relaxed";
end

savefile = outdir + algo + "_" + N + "_vars.mat";
save(savefile, 'xx', 'obj_originals', 'w_opt_vec', 'compl_vals', 'ind_x', 'ind_z', 'N', 'h', 'time_vals', 'rho_vals', 'nxsamples');

end