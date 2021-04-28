function [] = PrintStats(k, i, stat, compl, R, alphak, Mk, feasible, norm_pk, iterQPO)
logicalStr = {'false', ' true'};

% Print some stats
if (mod(i, 10) == 1)
    fprintf("------|-------|--------------|-------------|---------|-------------|-----------|----------|-----------|--------\n");
    fprintf("    k |     i | stationarity | complement. | penalty | step length | merit val | feasible |  norm pk  | qpIter \n");
    fprintf("------|-------|--------------|-------------|---------|-------------|-----------|----------|-----------|--------\n");
end
fprintf("%5d | %5d | %12.2g | %11.2g | %7.2g | %11.2g | %9.6g |    %s | %9.2g | %6d \n", k, i, stat, compl, R, alphak, Mk, logicalStr{feasible + 1}, norm_pk, iterQPO);
end