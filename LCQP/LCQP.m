function [ primalSteps, dualSteps, fVals, RVals, exitflag, auxOutput ] = LCQP(Q, g, A, lb, ub, lbA, ubA, L, R, params)
%% Is input consistent?
[nx, nc, ncomp] = PerformInputConsistencyTest(Q, g, A, lb, ub, lbA, ubA, L, R);

% Important for dual variables
A = [A; L; R];
lbA = [lbA; zeros(2*ncomp,1)];
ubA = [ubA; inf(2*ncomp,1)];

%% Add unbounded box constraints if none are passed
if (isempty(lb))
    lb = -inf(nx,1);
end
if (isempty(ub))
    ub = inf(nx,1);
end

%% Parameter handling
% qpOASES mode
options = qpOASES_options('default');
if (isfield(params, 'qpOASES_options'))
    options = params.qpOASES_options;
end

% Tolerence for convergence criterion
complementarityTolerance = 1e-10;
if (isfield(params, 'complementarityTolerance'))
    complementarityTolerance = params.complementarityTolerance;
end

% Tolerence for convergence criterion
stationarityTolerance = options.terminationTolerance*10;
if (isfield(params, 'stationarityTolerance'))
    stationarityTolerance = params.stationarityTolerance;
end

xk = zeros(size(nx,1));
% Initial starting point
if (isfield(params, 'x0'))
    % Initialize primals
    xk = params.x0;
else
    warning("No initial point given.");
end

useInitializationStruct = false;
% Initialization struct
if (isfield(params, 'initializationStruct'))
    initializationStruct = params.initializationStruct;
    useInitializationStruct = true;
    
    if (length(initializationStruct.guessedWorkingSetC) == nc)
        tmp = initializationStruct.guessedWorkingSetC;
        initializationStruct.guessedWorkingSetC = [tmp; zeros(2*ncomp,1)];
    end
end

solveZeroPenaltyFirst = true;
% Solve initial QP with or without penalization
if (isfield(params, 'solveZeroPenaltyFirst'))
    solveZeroPenaltyFirst = params.solveZeroPenaltyFirst;
end

% Store all steps?
storeSteps = false;
if (isfield(params, 'storeSteps'))
    storeSteps = params.storeSteps;
end

% Initial penalty parameter (check if given else 0)
Rk = 1;
if (isfield(params, 'R0'))
    Rk = params.R0;
end

% Update function for penalty parameter
penaltyUpdater = @(R) 2*R;
if (isfield(params, 'penaltyUpdater'))
    penaltyUpdater = params.penaltyUpdater;
end

Rbreak = 1000000;
% Penalty parameter value required to break the iteration.
if (isfield(params, 'Rbreak'))
    Rbreak = params.Rbreak;
end

maxIter = 100000;
% Penalty parameter value required to break the iteration.
if (isfield(params, 'maxIter'))
    maxIter = params.maxIter;
end

options.maxIter = maxIter;

printStats = false;
% Initial starting point
if (isfield(params, 'printStats'))
    printStats = params.printStats;
    assert(islogical(printStats));
end

useSparseMatrices = false;
% Initial starting point
if (isfield(params, 'useSparseMatrices'))
    useSparseMatrices = params.useSparseMatrices;
    assert(islogical(useSparseMatrices));
end

if (useSparseMatrices)
    Q = sparse(Q);
    A = sparse(A);
    L = sparse(L);
    R = sparse(R);
end

%% Helpers
% global iter counter
iter = 0;

% Update iteration whenever x is updated
k = 0;

% Originial objective
F = @(primal) 1/2*primal'*Q*primal + g'*primal;

% Symmetric penalty matrix
C = L'*R + R'*L;

% Penalty function 
phi = @(primal) 1/2*primal'*C*primal;
grad_phi = @(primal) C*primal;

% Stationarity and complementarity
Stationarity = @(primal, pen, dual_xk, dual_gk) Q*primal + g + pen*grad_phi(primal) - dual_xk - A'*dual_gk;

% Merit function
merit_fun = @(primal, pen) F(primal) + pen*phi(primal);

% Initial return values
primalSteps = xk;
RVals = Rk;
fVals = [];
dualSteps = [];

%% Initiate qpOASES homotopy
d_old = grad_phi(xk);

if (solveZeroPenaltyFirst)
    QP = qpOASES_sequence('i', Q, g, A, lb, ub, lbA, ubA, options, params.x0);
    [xk,~,exitflag,iterQPO,lk,auxOutput] = qpOASES_sequence('h', QP, g + Rk*d_old, lb, ub, lbA, ubA, options);     
else
    [QP,xk,~,exitflag,iterQPO,lk,auxOutput] = qpOASES_sequence('i', Q, g + Rk*d_old, A, lb, ub, lbA, ubA, options);
end
    
if (exitflag ~= 0)
    qpOASES_sequence('c', QP);
    return;
end

lam_xk = lk(1:length(xk));
lam_gk = lk((length(xk)+1):end);

alpha = 1;
pk = 1;

i = 0;
% Inner SQP loop
while ( true )        
    i = i + 1;

    iter = iter + 1;

    %% Store steps and update output
    if (storeSteps)
        primalSteps = [primalSteps, xk];
        dualSteps = [dualSteps, lk];
        fVals = [fVals, F(xk)];
        RVals = [RVals; Rk];
    else
        primalSteps = xk;
        dualSteps = lk;
        fVals = F(xk);
        RVals = Rk;
    end    


    %% Termination checks       
    % Evaluate stationarity and complementarity
    stat = norm(Stationarity(xk, Rk, lam_xk, lam_gk));

    % Iteration stats
    if (printStats)
        compl = abs(phi(xk));
        isFeasible = PrimalDualFeasibilitySequence(lb, ub, lbA, ubA, A*xk, xk, lam_xk, lam_gk, options.terminationTolerance);
        PrintStats(k, i, stat, compl, Rk, alpha, merit_fun(xk, Rk), isFeasible, norm(pk), iterQPO);
    end
    
    % Need this alternative stationarity check due to qpOASES termination
    % criterion
    if (stat <= options.terminationTolerance*50)
        if (Rk > Rbreak)
            fprintf("Exceeded maximum penalty value!\n"); 
            qpOASES_sequence('c', QP);
            return;
        end
        compl = abs(phi(xk));
        if (compl < complementarityTolerance)
            qpOASES_sequence('c', QP);       
            
            L_indices = nx+nc+1:nx+nc+ncomp;
            R_indices = nx+nc+ncomp+1:nx+nc+2*ncomp;
            lk(L_indices) = lk(L_indices) - Rk*R*xk;
            lk(R_indices) = lk(R_indices) - Rk*L*xk;         
            
            % GetStationarityType(xk, lk, L, R, L_indices, R_indices, complementarityTolerance);
            
            return;
        end
        % Update penalty parameter
        k = k + 1;
        i = 0;
        Rk = penaltyUpdater(Rk);
    end
    if (i > maxIter)
        disp('Maximum number of iterations reached. Exiting without convergence');
        
        exitflag = 2;
        qpOASES_sequence('c', QP);
        return;
    end
    
    %% Step computation
    d = grad_phi(xk) + (rand(nx, 1) - 0.5)*4*eps;    

    [x_new,~,exitflag,iterQPO,lambda,auxOutput] = qpOASES_sequence('h', QP, g + Rk*d, lb, ub, lbA, ubA, options);      
    
    if (exitflag ~= 0)
        qpOASES_sequence('c', QP);
        return;
    end

    %% Otpimal step length computation
    pk = x_new - xk;
    alpha = GetOptimalStepLength(C, g, Rk, Q, pk, xk);
    
    % Handle Optimal 0 step
    if (alpha == 0)
        if (norm(pk) < stationarityTolerance || phi(xk) < complementarityTolerance)
            alpha = 1;
        else
            % Update penalty parameter
            k = k + 1;
            i = 0;
            Rk = penaltyUpdater(Rk);
        end
    end
    
    %% Step update
    xk = xk + alpha*pk;
    lam_xk = lambda(1:nx);
    lam_gk = lambda((nx+1):end);
    lk = [lam_xk; lam_gk];       
end
end