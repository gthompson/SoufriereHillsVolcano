function plot_rsam_event_thresholds();
% Author: Glenn Thompson 2001
% Description: Plots RSAM event thresholds for each station
% Usage: plot_rsam_event_thresholds;

% Open new figure window
figure;

% Load RSAM alarm data file
load rsam_alarms.mat

% plot data
height=0.9/(length(stations));
styles={'b-';'k-';'g-.';'r-';'m-.';'c-'};
l=length(dnum_e);
for i=1:length(stations)
   	semilogy(dnum_e,threshold_e(:,i),styles{i});
   	set(gca,'XLim',[datenum(2000,1,1) datenum(2000,12,31)]);
   	set(gca,'YLim',[40 2000]);
    hold on;
end
title('RSAM alarm - event thresholds');
legend(stations,0);
hold off;
datetick('x',12);

