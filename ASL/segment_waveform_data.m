% break 1-day waveform files - which are too big for my laptop - into 10
% minute waveforms
clear all
ds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/analysis_events/data/%s.%s.%s.%s.D.%04d.%03d.miniseed','network','station','location','channel','year','jday')
ChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHE');
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHN')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHZ')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'10','HDF')];
snum(1) = datenum(2018,4,3,1,0,0); % Julian day 093
enum(1) = datenum(2018,4,3,4,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
snum(2) = datenum(2018,5,19,18,0,0); % Julian day 139
enum(2) = datenum(2018,5,19,20,0,0); % 17-20*3, 10-11
snum(3) = datenum(2018,9,25,9,0,0); % Julian day 268
enum(3) = datenum(2018,9,25,11,30,0); % 12-20, 8-11*2
snum(4) = datenum(2018,10,15,9,30,0); % Julian day 288 
enum(4) = datenum(2018,10,15,11,30,0); % 2-3, 10-12
snum(5) = datenum(2018,11,12,14,0,0); % Julian day 316  % distromidor sees rain up to 40mm/hr
enum(5) = datenum(2018,11,12,19,0,0); % 13-24
% sandmining?
snum(6) = datenum(2018,4,3,13,0,0); % Julian day 093
enum(6) = datenum(2018,4,3,16,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
snum(7) = datenum(2018,4,3,17,0,0); % Julian day 093
enum(7) = datenum(2018,4,3,20,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
iceweb.usf_calibrations;
for eventnum=1:numel(snum)
    for ctagnum=1:numel(ChannelTagList)
        ctag = ChannelTagList(ctagnum);
        w = waveform(ds, ctag, snum(eventnum), enum(eventnum));
        wavfilename = sprintf('%s_%s.mat',ctag.string(),datestr(snum(eventnum),30) );
        disp(wavfilename);
        save(wavfilename,'w');
    end
end
% % 
% % defaultSamplingIntervalSeconds = 60;
% % products.waveform_plot.doit = true;
% % products.rsam.doit = false;
% % products.rsam.samplingIntervalSeconds = [1 defaultSamplingIntervalSeconds]; % [1 60] means record RSAM data at 1-second and 60-second intervals
% % products.rsam.measures = {'max';'mean';'median'}; % {'max';'mean';'median'} records the max, mean and median in each 1-second and 60-second interval
% % products.spectrograms.doit = false; % whether to plot & save spectrograms
% % products.spectrograms.plot_metrics = true; % superimpose metrics on spectrograms
% % products.spectrograms.timeWindowMinutes = 10; % 60 minute spectrograms. 10 minute spectrograms is another common choice
% % products.spectrograms.fmin = 0.5;
% % products.spectrograms.fmax = 100; % Hz
% % products.spectrograms.dBmin = 60; % white level
% % products.spectrograms.dBmax = 120; % pink level
% % products.spectral_data.doit = true; % whether to compute & save spectral data
% % products.spectral_data.samplingIntervalSeconds = defaultSamplingIntervalSeconds; % DO NOT CHANGE! spectral data are archived at this interval
% % products.soundfiles.doit = false;
% % products.helicorders.doit = false;
% % products.reduced.doit = false;
% % products.reduced.samplingIntervalSeconds = defaultSamplingIntervalSeconds;
% % products.removeWaveformFiles = false;
% % 
% % daily plots
% % products.daily.spectrograms = false;
% % products.daily.helicorders = false;
% % products.daily.rsamplots = false;
% % products.daily.spectralplots = false;
% % 
% % generate list of timewindows
% % nummins = 10;
% % timewindows = iceweb.get_timewindow(enum, nummins, snum);
% % 
% % loop over timewindows
% % for count = 1:length(timewindows.start)
% %     hh = datestr(timewindows.start(count),'HH');
% %     mm = datestr(timewindows.start(count),'MM');
% %     if strcmp(hh,'00') && strcmp(mm, '00') || count==1
% %         fprintf('\n%s ',datestr(timewindows.start(count),26));
% %     end
% %     if strcmp(mm,'00')
% %         fprintf('%s ',hh);
% %     end
% %     iceweb.process_timewindow(products_dir, networkName, ...
% %         subnetName, ChannelTagList, timewindows.start(count), ...
% %         timewindows.stop(count), ds, products, calibObjects);
% % 
% % end