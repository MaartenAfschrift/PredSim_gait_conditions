function [] = SetDescriptiveStatistics_Fig(h, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

regmethod = 'standard';
if ~isempty(varargin)
    regmethod= varargin{1};
end

% get descriptive statistics:
all_axes = findobj(h, 'Type', 'axes');
for i =1:length(all_axes)
    target_subplot = all_axes(i);
    lines = findobj(target_subplot, 'Type', 'line');
    xdat = [];
    ydat = [];
    for j = 1:length(lines)
        xdat = [xdat; lines(j).XData'];
        ydat = [ydat; lines(j).YData'];
    end
    % delete min and max values (lines)
    iDel = find(xdat == min(xdat) | xdat == max(xdat) | xdat == 0 | isnan(xdat));
    xdat(iDel) = [];
    ydat(iDel) = [];

    % compute correlation
    if strcmp(regmethod,'standard')
        stats = regstats(xdat,ydat,'linear');
        text(target_subplot,nanmean(xdat)+nanstd(xdat),nanmean(ydat), {['Rsq ' num2str(round(stats.rsquare,2))];...
            ['coeff ' num2str(round(stats.beta(2),2))]; ['rmse ' num2str(round(rmse(stats.yhat,ydat),2))]});
    elseif strcmp(regmethod,'paper_predint')
        % stats for correlation coefficient
        stats = regstats(xdat,ydat,'linear');
        % coefficients based on zero intercept (as we compute delta values)
        coeff = xdat\ydat;
        % rmse is difference between xdat and ydat
        rmse_out = rmse(xdat,ydat);
        text(target_subplot,nanmean(xdat)+nanstd(xdat),nanmean(ydat), {['Rsq ' num2str(round(stats.rsquare,2))];...
            ['coeff ' num2str(round(coeff,2))]; ['rmse ' num2str(round(rmse_out,2))]});

    end
end


end