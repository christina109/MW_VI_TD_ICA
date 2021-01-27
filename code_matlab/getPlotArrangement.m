function [nrow, ncol] = getPlotArrangement(nPlot)
    
if nPlot <= 4
    nrow = 1;
    ncol = nPlot;
elseif nPlot <= 9
    ncol = 3;
    nrow = ceil(nPlot/3);
elseif nPlot <= 16
    ncol = 4;
    nrow = ceil(nPlot/4);
elseif nPlot <= 25
    ncol = 5;
    nrow = ceil(nPlot/5);
elseif nPlot <= 36
    ncol = 6;
    nrow = ceil(nPlot/6);
elseif nPlot <= 49
    ncol = 7;
    nrow = ceil(nPlot/7);
elseif nPlot <= 64
    ncol = 8;
    nrow = ceil(nPlot/8);
end

end
