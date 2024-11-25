function [l] = PlotIndexPlot(IndexPlot,Table,Headers,Cs,mk,xHeader,yHeader,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% plot difference values
PlotDelta = true;
if ~isempty(varargin)
    PlotDelta = varargin{1};
end


for i=1:length(IndexPlot)
    iSel = IndexPlot(i).iSel; iRef = IndexPlot(i).iRef;
    if PlotDelta
        DeltaDat = Table(iSel,:) - Table(iRef,:);
        l = plot(DeltaDat(:,strcmp(Headers,xHeader)),...
            DeltaDat(:,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);
    else
        l = plot(Table(iSel,strcmp(Headers,xHeader)),...
            Table(iSel,strcmp(Headers,yHeader)),...
            'o','Color',Cs,'MarkerFaceColor',Cs,'MarkerSize',mk);

    end
end
end