function rsam_alarms_menu();
% Author: Glenn Thompson
% Description: A menu/wrapper for RSAM alarm functions
% Usage: rsam_alarms;

warning off;
choice = 0;
while choice<7
    choice = menu('RSAM ALARMS','Update Matlab RSAM alarm datafile','Export RSAM alarm data to Excel','Plot event thresholds','Text display RSAM alarms last week','Histogram of weekly alarms','Close all figures','Exit');
    switch choice
        case 1, update_rsam_alarms;
        case 2, rsam_alarm2excel;
        case 3, plot_rsam_event_thresholds;
        case 4, text_rsam_alarms(now-7-4/24,now-4/24);
        case 5, histogram_rsam_alarms();
        case 6, close all;
    end
end