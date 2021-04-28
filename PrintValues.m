function [] = PrintValues(plotStructs, names, key)

n = 0; m = 0;

for i=1:length(plotStructs)
    plotStruct = plotStructs{i};
    [tn, tm] = size(plotStruct(key));
    n = max(n, tn);
    m = max(m, tm);
end

fprintf("     method      |    " + key + "\n");
fprintf("-----------------|-----------------\n");

for i=1:length(plotStructs)
    plotStruct = plotStructs{i};
    
    if (all(size(plotStruct(key)) == [n m]))
        fprintf("%16s | %12.2g\n", names{i}, mean(reshape(abs(plotStruct(key)), n*m, 1)));
    end
end

fprintf("-----------------|-----------------\n\n");

end