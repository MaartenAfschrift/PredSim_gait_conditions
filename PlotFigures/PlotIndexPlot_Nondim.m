function [l] = PlotIndexPlot_Nondim(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,scale_exp,scale_sim,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% plot difference values
PlotDelta = true;
if ~isempty(varargin)
    PlotDelta = varargin{1};
end

% scale variable with specific number
scale = 1;
if length(varargin)>1
    scale = varargin{2};
end


for i=1:length(IndexPlot)
    iSel = IndexPlot(i).iSel; 
    iRef = IndexPlot(i).iRef;
    if PlotDelta
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        xdat = DeltaDat(:,strcmp(Headers,xHeader));
        ydat = DeltaDat(:,strcmp(Headers,yHeader));
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*scale,ydat*scale,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    else
        xdat = Table(iSel,strcmp(Headers,xHeader));
        ydat = Table(iSel,strcmp(Headers,yHeader));
        xdat = xdat./scale_exp.*scale_sim;
        l = plot(xdat*scale,ydat*scale,'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

    end
end
end