close all; clear all; clc;

%% Add required paths
% Adjust CasADi and qpOASES path to your local settings:
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
mpcc_modes = [2 3 4 5];

% In order to compare the sparse schur method against the dense method with
% sparse matrices it is recommended to adjust the experiment name below to
% <expname>/schur. This behavior is supported in the creation of
% performance plots.

% Number of initializations tested
nxsamples = 100;

% Number of discretization nodes (scales the IVOCP)
N_vals = 50:5:150;

% Experiment name and output directory
expname = "experiment";
outdir = "saved_variables/" + expname + "/";

% Create output directory if it doesn't exist
if ~exist(outdir, 'dir')
   mkdir(outdir)
end

% Penalty homotopy or regularization homotopy
rho_vals = [0 0.1 0.2 0.4 0.8 1.6];
sigma_vals = [1e-1 5e-2 1e-2 5e-3 1e-3];

for i = 1:length(N_vals)
    for j = 1:length(mpcc_modes)
        fprintf("Running [mode = %d, N = %d]\n", mpcc_modes(j), N_vals(i));
        
        if (mpcc_modes(j) == 4 || mpcc_modes(j) == 5)
            % Regularization
            IVOCP(mpcc_modes(j), sigma_vals, nxsamples, N_vals(i), outdir);
        else
            % Penaltization
            IVOCP(mpcc_modes(j), rho_vals, nxsamples, N_vals(i), outdir);                
        end
    end
end

%% Generate output
close all; clear all; clc;
CreatePerformancePlots(outdir);
