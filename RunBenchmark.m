%% Add required paths
% Adjust CasADi and qpOASES path!
addpath("~/Documents/MATLAB/casadi-linux-matlabR2014b-v3.5.5");
addpath("~/Documents/qpOASES/3.2/interfaces/matlab");

% Shipped directories
addpath("LCQP");
addpath("Collocation");

%% Settings
% This script can be used to run the solve the IVOCP with the different 
% solution variants (mpcc_mode).
% 1: LCQP
% 2: LCQP Sparse
% 3: IPOPT Penalty
% 4: IPOPT Smoothed
% 5: IPOPT Relaxed
mpcc_modes = 2; %[2 3 4 5];

% Number of initializations tested
nxsamples = 100;

% Number of discretization nodes (scales the IVOCP)
N_vals = 50:155:150;

% Experiment name and output directory
expname = "test";
outdir = "saved_variables/" + expname;

% Create output directory if it doesn't exist
if ~exist(outdir, 'dir')
   mkdir(outdir)
end

% Penalty homotopy or regularization homotopy
rho_vals = [0 0.1 0.2 0.4 0.8 1.6];
sigma_vals = [1e-1 5e-2 1e-2 5e-3 1e-3];

for i = 1:length(N_vals)
    for j = 1:length(mpcc_modes)
        if (mpcc_modes(j) == 4 || mpcc_modes(j) == 5)
            % Regularization
            IVOCP(mpcc_modes(j), sigma_vals, nxsamples, N_vals(i), outdir);
        else
            % Penaltization
            IVOCP(mpcc_modes(j), rho_vals, nxsamples, N_vals(i), outdir);                
        end

        fprintf("i = %d | j = %d\n", i, j);
    end
end

% Generate output
IVOCP_Performance_Plots(outdir);
