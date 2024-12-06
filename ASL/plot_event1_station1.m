clear all
close all
snum(1) = datenum(2018,4,3,1,0,0); % Julian day 093
enum(1) = datenum(2018,4,3,4,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
snum(2) = datenum(2018,5,19,18,0,0); % Julian day 139
enum(2) = datenum(2018,5,19,20,0,0); % 17-20*3, 10-11
snum(3) = datenum(2018,9,25,9,0,0); % Julian day 268
enum(3) = datenum(2018,9,25,11,30,0); % 12-20, 8-11*2
snum(4) = datenum(2018,10,15,9,30,0); % Julian day 288 
enum(4) = datenum(2018,10,15,11,30,0); % 2-3, 10-12
snum(5) = datenum(2018,11,12,14,0,0); % Julian day 316  % disdrometer sees rain up to 40mm/hr
enum(5) = datenum(2018,11,12,19,0,0); % 13-24
% sandmining?
snum(6) = datenum(2018,4,3,13,0,0); % Julian day 093
enum(6) = datenum(2018,4,3,16,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
snum(7) = datenum(2018,4,3,17,0,0); % Julian day 093
enum(7) = datenum(2018,4,3,20,0,0); % 0-4*4, 10-11, 12-20 (14-15, 17-18)
ChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHE');
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHN')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHZ')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'10','HDF')];
rsambinsize = 10; % seconds
%% each event, by station
close all
stations = unique({ChannelTagList.station});
for eventnum=1:numel(snum)
    for stanum = 1:numel(stations)
        thissta = stations{stanum};
        ctags = ChannelTagList(matching(ChannelTagList, sprintf('MV.%s..',thissta)))
        wvector = [];
        for ctagnum=1:numel(ctags)
            ctag = ctags(ctagnum);
            wavfilename = sprintf('%s_%s.mat',ctag.string(),datestr(snum(eventnum),30) );
            disp(wavfilename);
            load(wavfilename);
            wvector = [wvector w];
        end
        calibObjects = iceweb.usf_calibrations(ctags);
        wvector = iceweb.apply_calib(wvector, calibObjects)
        plot_panels(wvector);
        print('-dpng',sprintf('event%d_%s.png',eventnum,ctag.station));
        
%         % integrate to displacement / pressure
%         filterfreqs = [0.01 10];
%         wdetrend = detrend(wvector);
%         wtaper = taper(wdetrend,0.1);
%         wclean = clean(wtaper, filterfreqs);
%         
%         
%         plot_panels(wclean);
%         
%         wdisp = wclean;
%         for c=1:numel(wclean)
%             mychan = get(wclean(c),'channel');
%             if mychan(2) == 'H' % integrate
%                 wdisp(c) = integrate(wclean(c),'trapz');
%             end
%         end
%         plot_panels(wdisp);
%         print('-dpng',sprintf('event%d_%s_displacement.png',eventnum,ctag.station));
%         %anykey = input('any key')
        
        r = waveform2rsam(wvector, 'mean', rsambinsize);
        r.plot2();
        print('-dpng',sprintf('event%d_%s_rsam_mean.png',eventnum,ctag.station));


        r = waveform2rsam(wvector, 'median', rsambinsize);
        r.plot2();    
        print('-dpng',sprintf('event%d_%s_rsam_median.png',eventnum,ctag.station));

        
        
        close all
    end
    
end

%% each event, by channel
close all
channels = unique({ChannelTagList.channel});
for eventnum=1:numel(snum)
    for channum = 1:numel(channels)
        thischan = channels{channum};
        ctags = ChannelTagList(matching(ChannelTagList, sprintf('MV...%s',thischan)))
        wvector = [];
        for ctagnum=1:numel(ctags)
            ctag = ctags(ctagnum);
            wavfilename = sprintf('%s_%s.mat',ctag.string(),datestr(snum(eventnum),30) );
            disp(wavfilename);
            load(wavfilename);
            wvector = [wvector w];
        end
        
        % true ground velocity / pressure
        calibObjects = iceweb.usf_calibrations(ctags);
        wvector = iceweb.apply_calib(wvector, calibObjects);
        plot_panels(wvector);
        print('-dpng',sprintf('event%d_%s.png',eventnum,ctag.channel));
        
%         % integrate to displacement / pressure
%         wclean = clean(wvector);
%         wdisp = wclean;
%         for c=1:numel(wclean)
%             mychan = get(wclean(c),'channel');
%             if mychan(2) == 'H' % integrate
%                 wdisp(c) = integrate(wclean(c));
%             end
%         end
%         plot_panels(wdisp);
%         print('-dpng',sprintf('event%d_%s_displacement.png',eventnum,ctag.channel));
%         %anykey = input('any key')

        r = waveform2rsam(wvector, 'mean', rsambinsize);
        r.plot2();
        print('-dpng',sprintf('event%d_%s_rsam_mean.png',eventnum,ctag.channel));


        r = waveform2rsam(wvector, 'median', rsambinsize);
        r.plot2();    
        print('-dpng',sprintf('event%d_%s_rsam_median.png',eventnum,ctag.channel));



        close all
    end
end

%% All events for each chantag for  events
close all
clear w r wvector ctags
channels = unique({ChannelTagList.channel});

for ctagnum = 1:numel(ChannelTagList)
    ctag = ChannelTagList(ctagnum);
    wvector =[];
    for eventnum=1:5
        wavfilename = sprintf('%s_%s.mat',ctag.string(),datestr(snum(eventnum),30) );
        disp(wavfilename);
        load(wavfilename);
        wvector = [wvector w];
    end

    % true ground velocity / pressure
    calibObjects = iceweb.usf_calibrations(ctag);
    wvector = iceweb.apply_calib(wvector, calibObjects);
    plot_panels(wvector,'alignWaveforms',1);
    print('-dpng',sprintf('event%d_%s.png',eventnum,ctag.string()));

%     % integrate to displacement / pressure
%     wclean = clean(wvector);
%     wdisp = wclean;
%     for c=1:numel(wclean)
%         mychan = get(wclean(c),'channel');
%         if mychan(2) == 'H' % integrate
%             wdisp(c) = integrate(wclean(c));
%         end
%     end
%     plot_panels(wdisp);
%     print('-dpng',sprintf('event%d_%s_displacement.png',eventnum,ctag.string()));
%     %anykey = input('any key')

    r = waveform2rsam(wvector, 'mean', rsambinsize);
    r.plot2();
    print('-dpng',sprintf('event%d_%s_rsam_mean.png',eventnum,ctag.string()));

    
    r = waveform2rsam(wvector, 'median', rsambinsize);
    r.plot2();    
    print('-dpng',sprintf('event%d_%s_rsam_median.png',eventnum,ctag.string()));


    close all
end

% next stage should be to make 1-s and 60-s RSAM and then export
% as simple datenum and amplitude structs for Dav to plot

