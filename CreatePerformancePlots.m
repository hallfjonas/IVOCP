function[] = CreatePerformancePlots( outdir )
%% Load variables
lcqp_sparse_mat_Files = dir(outdir + "LCQP_Sparse_*_vars.mat");
schur_lcqp_sparse_mat_Files = dir(outdir + "schur/LCQP_Sparse_*_vars.mat");
nlp_pen_Files = dir(outdir + "nlp_pen_*_vars.mat");
relaxed_Files = dir(outdir + "relaxed_*_vars.mat");
smoothed_Files = dir(outdir + "smoothed_*_vars.mat");

% Variable names
varNames = [...
    "xx", "obj_originals", "w_opt_vec", "compl_vals", "ind_x", ...
    "ind_z", "N", "h", "time_vals", "rho_vals", "nxsamples", 
];

% create container for each method
lcqpSparse = containers.Map;
schurLcqpSparse = containers.Map;
nlp_pen = containers.Map;
smoothed = containers.Map;
relaxed = containers.Map;
lcqpSparse("files") = lcqp_sparse_mat_Files;
schurLcqpSparse("files") = schur_lcqp_sparse_mat_Files;
nlp_pen("files") = nlp_pen_Files;
smoothed("files") = smoothed_Files;
relaxed("files") = relaxed_Files;

% Store all containers in one container
methods = containers.Map;
methods("LCQP") = lcqpSparse;
methods("LCQP Schur") = schurLcqpSparse;
methods("IPOPT Pen") = nlp_pen;
methods("Relaxed") = relaxed;
methods("Smoothed") = smoothed;

% Fill containers
mKeys = methods.keys;
for k = 1:methods.Count
    c = methods(mKeys{k});
    for i = 1:length(varNames)
        c(varNames(i)) = [];
    end
    
    allFiles = c("files");

    for i=1:length(allFiles)
        load(allFiles(i).folder + "/" + allFiles(i).name)
        
        w_opt_vec = w_opt_vec(2,:)';
        
        for j = 1:length(varNames)
            c(varNames(j)) = [c(varNames(j)) eval(varNames(j))];
        end
    end
    
    % Sort by number of variables
    [~, IND] = sort(c("N"));
    keys = c.keys;
    for i = 1:c.Count
        tmp = c(keys{i});
        if (size(tmp,2) == 1)
            c(keys{i}) = tmp(IND);
        else
            c(keys{i}) = tmp(:,IND);
        end
    end
end

% Create sorted plot structure containing all methods
methodLegendEntries = {...
    'LCQP', 'LCQP Schur', ...
    'IPOPT Pen', 'IPOPT Smoothed', 'IPOPT Relaxed'...
};
plotStructs = { ...
    methods("LCQP"), ...
    methods("LCQP Schur"), ...
    methods("IPOPT Pen"), ...
    methods("Smoothed"), ...
    methods("Relaxed") ...
};

%% Plot settings
% Globally use latex interpreter
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Colors for generic lengths
% colors = distinguishable_colors(length(methods));
% markers = {'d','s','o','<','x','v','+','*','^','>'};
% linetypes = {'-', '--', '-', '--','-.', ':', ':'};

% Fixed chosen colors and markers
%colors = ['k', 'r', 'b', 'm', 'y', 'g'];
colors = distinguishable_colors(5);
markers = {'d','x','s','>', 'o', '+'};
linetypes = {'-', '--', ':', '-.','-.', '--'};

%% Compute Analytic Solution for comparison
x = linspace(-1.9,-0.9,1000);
L = zeros(length(x), 1);
T = 2;
for i = 1:length(x)
    x0 = x(i);
    ts = -x0/3;
    
    c0 = (T+(x0-5)/3)^2;
    c1 = 8/3*ts^3 +8/3*x0*ts^2 + 8/9*x0^2*ts +1/3*T^3+1/3*x0*T^2+x0^2*T/9;
    L(i) = c0 + c1;    
end

% Get optimal x0*
[~, ind_min] = min(L);
x0_opt = x(ind_min);

%% Figure 1: Quality Plots
figure(1); box on; hold on;
for  i=1:length(plotStructs)
    plotStruct = plotStructs{i};
    semilogy(...
        plotStruct("N"), mean(abs(plotStruct("w_opt_vec")-x0_opt), 1), ...
        'Marker', markers{i}, 'LineStyle', linetypes{i}, ...
        'Color', colors(i,:), "DisplayName", methodLegendEntries{i} ...
    );
end
grid on;

% ylim([0 0.55]);
legend("Location", "northwest");
xlabel("Number of discretization nodes");
ylabel("$\mathrm{mean}(|x_0^\ast - \hat{x}_0|)$ per $100$ OCPs");

% Styling
lines = findobj(gcf,"Type","Line");
for i = 1:numel(lines)
    lines(i).LineWidth = 1.3;
end

% Export using Matlab2Tikz and PNG
filename = 'QualityPlot';
relativeImagesPath = outdir;
saveas(gcf,strcat(relativeImagesPath, filename, '.png'));

% Uncomment if exporting to tikz if desired (must be installed)
% matlab2tikz(strcat(relativeImagesPath, filename, '.tikz'));

%% Figure 2: Time plots
figure(2); box on; hold on;
for  i=1:length(plotStructs)
    plotStruct = plotStructs{i};
    semilogy(...
        plotStruct("N"), mean(plotStruct("time_vals"), 1), ...
        'Marker', markers{i}, 'LineStyle', linetypes{i}, ...
        'Color', colors(i,:), "DisplayName", methodLegendEntries{i}...
    );
end


% Styling
legend("Location", "southeast");
xlabel("Number of discretization nodes");
ylabel("Average CPU time $[s]$");
grid on;
lines = findobj(gcf,"Type","Line");
for i = 1:numel(lines)
    lines(i).LineWidth = 1.3;
end

% Set to log scale
set(gca,'yscale','log')
 
% Export using Matlab2Tikz and PNG
filename = 'TimePlot';
relativeImagesPath = outdir;
saveas(gcf,strcat(relativeImagesPath, filename, '.png'));

% Uncomment if exporting to tikz if desired (must be installed)
% matlab2tikz('TimePlot.tikz');

%% Complementarity and Quality Averages
for i = 1:length(plotStructs)
    plotStruct = plotStructs{i};
    plotStruct("x0_diff") = abs(plotStruct("w_opt_vec") - x0_opt);
    plotStructs{i} = plotStruct;
end

PrintValues(plotStructs, methodLegendEntries, "x0_diff");
PrintValues(plotStructs, methodLegendEntries, "compl_vals");

end
