% These are the lahars spotted by Paddy. Need to:
% 1. Get the data from FTP server
%    - done 2019/03/24
% 2. Get the iceweb products from bigPC (or newton)
%    - done 2019/03/24
% 3. Load waveform data into GISMO
% 4. Plot panels for all seismic and all infrasound
% 5. Compute 1-s RSAM
% 6. Try to locate
%datasourceObject = datasource('miniseed', '/media/sdb1/belhamstudy/download2/%s/%04d/%s/%s/%s.D/%s.%s.%s.%s.D.%04d.%03d.miniseed','station','year','network','station','channel','network','station','location','channel','year','jday')
datasourceObject = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/data/%s.%s.%s.%s.D.%04d.%03d.miniseed','network','station','location','channel','year','jday')
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
PRODUCTS_TOP_DIR = './icewebproducts';
subnetName = 'Belham'; % can be anything
% output directory is PRODUCTS_TOP_DIR/network_name/subnetName
% get parameters. these control what iceweb does.
products.DAILYONLY = false; % set true to only produce day plots
iceweb.get_params;
products.soundfiles.doit = false;
products.helicorders.doit = false;
close all
% % for eventnum=1:numel(snum)
% %     disp(sprintf('Loading waveforms for event %d',eventnum));
% %     w = waveform(datasourceObject, ChannelTagList, snum(eventnum), enum(eventnum));
% %     w = iceweb.apply_calib(w, calibObjects);
% %     r = waveform2rsam(wtmp, 'median', 1.0);
% %     save(sprintf('laharevent%d_waveform.mat',eventnum));
% % end
for eventnum=1:numel(snum)
    %% run IceWeb  ********************** this is where stuff really happens
    iceweb.run_iceweb(PRODUCTS_TOP_DIR, subnetName, datasourceObject, ...
        ChannelTagList, snum(eventnum), enum(eventnum), gulpMinutes, products, calibObjects);
end



% %%
%     %
%     plot_panels(wtmp);
%     print(sprintf('laharevent%d_panelplot.png',eventnum),'-dpng');
%     %
%     s{eventnum} = waveform2rsam(wtmp, 'median', 1.0);
%     s{eventnum}.plot()
%     print(sprintf('laharevent%d_rsam.png',eventnum),'-dpng');
%     %
%     spectrogram(wtmp)
%     print(sprintf('laharevent%d_sgram.png',eventnum),'-dpng');
%     %
%     w = {w wtmp};
%     clear wtmp
%     save lahar_waveforms.mat
% end
%% rsam analysis of each lahar sequence
clc
for eventnum=1:5
    chans = {'HHZ';'HDF'};
    for channum = 1:numel(chans)
        chan = chans{channum};
        if chan(2)=='H'
            ctag = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHZ');
        elseif chan(2)=='D'
            ctag = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'10','HDF');
        end
        close all
        
        % get and plot rsam data
        filepattern = fullfile('/Volumes/data/Montserrat/fromBigPC/MV/Monty','SSSS.CCC.YYYY.MMMM.001.bob');
        r = iceweb.daily_rsam_plot(filepattern, floor(snum(eventnum)), ceil(enum(eventnum)), ctag, 'median');
        a=axis;
        set(gca,'YLim',[0 a(4)]);
        %legend(ctag.string());
        print( sprintf('%s_RSAMevent%d%s.png',datestr(snum(eventnum),'yyyymmdd'),eventnum,chan), '-dpng');
        
%         % existence plot
%         N = numel(r);
%         for counter=1:N
%             r2(counter)=r(counter);
%             data = r2(counter).data;
%             data(~isnan(data))=1+counter/100;
%             data(isnan(data))=0;
%             %data(2:2:end)=data(2:2:end)*(-1.0-counter/100);
%             r2(counter).data = data;
%         end
%         %r2.plot_panels()
%         r2.plot2();
%         legend(ChannelTagList.string());
%         a=axis;
%         set(gca,'YLim',[0 a(4)]);
%         pngfile = sprintf('%s_Existence_event%d%s.png',datestr(snum(eventnum),'yyyymmdd'),eventnum,chan);
%         print('-dpng',pngfile);
    end
end
%% rsam analysis
clc
startdate = datenum(2018,3,29);
enddate = datenum(2019,1,10); % last download
%enddate = datenum(2018,4,13); % testing
numweeks = ceil((enddate - startdate)/7);
for weeknum = 1 : numweeks
    chans = {'HHZ';'HDF'};
    for channum = 1:numel(chans)
        chan = chans{channum};
        if chan(2)=='H'
            ctag = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHZ');
        elseif chan(2)=='D'
            ctag = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'10','HDF');
        end
        snum = startdate + (weeknum-1)*7;
        enum = min([snum+7  enddate]);
        close all
        
        % get and plot rsam data
        filepattern = fullfile('/Volumes/data/Montserrat/fromBigPC/MV/Monty','SSSS.CCC.YYYY.MMMM.060.bob');
        r = iceweb.daily_rsam_plot(filepattern, snum, enum, ctag, 'max');
        a=axis;
        set(gca,'YLim',[0 a(4)]);
        legend(ChannelTagList.string());
        print( sprintf('RSAMplot_%s_week%03d.png',chan,weeknum), '-dpng');
        
        % existence plot
        N = numel(r);
        for counter=1:N
            r2(counter)=r(counter);
            data = r2(counter).data;
            data(~isnan(data))=1+counter/100;
            data(isnan(data))=0;
            %data(2:2:end)=data(2:2:end)*(-1.0-counter/100);
            r2(counter).data = data;
        end
        %r2.plot_panels()
        r2.plot2();
        legend(ChannelTagList.string());
        a=axis;
        set(gca,'YLim',[0 a(4)]);
        pngfile = sprintf('Existence_%s_week%03d.png',chan,weeknum);
        print('-dpng',pngfile);
    end
end



%%
% soh analysis
% 1. re-plot the RSAM data
% 2. Load soh data into MATLAB
% 3. Analyze
% 
% MTB2 data gaps from 03-Nov onwards
% MTB3 bad data from 22-July


SOHChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'','SOH');
SOHds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/soh/%s/%s.%s.D0.SOH.L.%04d.%03d.miniseed','station','network','station','year','jday')
%startdate = datenum(2018,3,29);
startdate = datenum(2019,3,19);
enddate = datenum(2019,3,26); % last download
numweeks = ceil((enddate - startdate)/7);
for weeknum = 1 : numweeks
    
    for stanum = 1:numel(SOHChannelTagList)
        close all
        ctag=SOHChannelTagList(stanum);
        snum = startdate + (weeknum-1)*7;
        if snum>=enddate
            break;
        end
        enum = min([snum+7  enddate]);
        wsoh = [];
        thisday = snum;
        while thisday< enum
            thisw = waveform(SOHds, ctag, thisday, thisday+1-1/86400)
            if ~isempty(wsoh)
                wsoh = combine([wsoh thisw]);
            else 
                wsoh = thisw;
            end
            thisday = thisday + 1;
        end

        if ~isempty(wsoh)
            figure
            set(gcf,'Position',[131 36 1236 772]);
            for c=1:6
                ax = subplot(3,2,c);
                h = plot(wsoh(c),'xunit','date','axeshandle', ax);
            end
            print( sprintf('sohplot1_%s_week%03d.png',ctag.station,weeknum), '-dpng');
            %
            figure
            set(gcf,'Position',[131 36 1236 772]);
            for c=7:9
                ax = subplot(3,1,c-6);
                h = plot(wsoh(c),'xunit','date','axeshandle', ax);
            end
            print( sprintf('sohplot2_%s_week%03d.png',ctag.station,weeknum), '-dpng');
            %
            figure
            set(gcf,'Position',[131 36 1236 772]);
            for c=10:15
                ax = subplot(3,2,c-9);
                h = plot(wsoh(c),'xunit','date','axeshandle', ax);
            end
            print( sprintf('sohplot3_%s_week%03d.png',ctag.station,weeknum), '-dpng');
        end
    end
end

% LCQ = clock quality. 100% = < 5us. 90% = <100us. 70%% < 200us (5000 Hz).


%%
% soh read files into MAT files
SOHChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'','SOH');
SOHds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/soh/%s/%s.%s.D0.SOH.L.%04d.%03d.miniseed','station','network','station','year','jday')
startdate = datenum(2018,3,29);
enddate = datenum(2019,3,26); % last download

close all
for stanum = 1:numel(SOHChannelTagList)
    ctag=SOHChannelTagList(stanum);
    wsoh = [];
    for snum=startdate:enddate
        try
            thisw = waveform(SOHds, ctag, snum, snum+1-1/86400);
            if ~isempty(wsoh)
                wsoh = combine([wsoh thisw]);
            else 
                wsoh = thisw;
            end
            wsoh
        save(sprintf('%s_%s_%s_soh.mat',ctag.station,datestr(startdate,'yyyymmdd'),datestr(snum,'yyyymmdd')),'wsoh')
        end
    end
end


%%

% 01 GAN = GPS antenna status 0=ok
% 02 GEL = GPS elevation
% 03 GLA = GPS lat
% 04 GLO = GPS long
% 05 GNS = # satellites
% 06 GPL = GPS PLL status (0=no lock)
% 07 GST = GPS status (0=off, 1=unlocked, 2=locked) 
% 08 LCE = difference between GPS and digitizer clock in us
% 09 LCQ = clock quality. 100% = < 5us. 90% = <100us. 70%% < 200us (5000 Hz).
% 10 VCO = VCO control voltage
% 11 VDT = digitizer temperature (millidegrees C)
% 12 VEC = digitizer current (mA)
% 13 VEI = input voltage (mV)
% 14 VM1 = mass position (uV)
% 15 VPB = digitizer buffer fullness (usually 100%)
close all
figure
for c=1:4
    station = sprintf('MTB%d',c);
    fname = sprintf('%s_20180329_20190326_soh.mat',station)
    load(fname);
    %get(wsoh,'channel')
    ax=subplot(4,1,c);
    hp=plot(wsoh(9),'xunit','date','axeshandle',ax)
    set(hp,'LineStyle','none','Marker','.')
    ylabel('clock quality %%');
    datetick('x')
end
%%
figure
for c=1:4
    station = sprintf('MTB%d',c);
    fname = sprintf('%s_20180329_20190326_soh.mat',station)
    load(fname);
    %get(wsoh,'channel')
    ax=subplot(4,1,c);
    hp=plot(wsoh(14)/1e6,'xunit','date','axeshandle',ax)
     set(hp,'LineStyle','none','Marker','.')
    ylabel('mass position V%');
    datetick('x')
end
% 
% 
%%
figure
for c=1:4
    station = sprintf('MTB%d',c);
    fname = sprintf('%s_20180329_20190326_soh.mat',station)
    load(fname);
    %get(wsoh,'channel')
    ax=subplot(4,1,c);
    hp=plot(wsoh(5),'xunit','date','axeshandle',ax)
    set(hp,'LineStyle','none','Marker','.')
    ylabel(sprintf('# Satellites'));
    datetick('x')
end
%
%%
figure
for c=2
    station = sprintf('MTB%d',c);
    fname = sprintf('%s_20180329_20190326_soh.mat',station)
    load(fname);
    wsoh(13)
    continue;
    %get(wsoh,'channel')
    ax=subplot(4,1,c);
    hp=plot(wsoh(13)/1e3,'xunit','date','axeshandle',ax)
    set(hp,'LineStyle','none','Marker','.')
    set(ax,'YLim',[9 15]);
    ylabel(sprintf('Voltage (V)'));
    datetick('x')
end

%%
% soh read files and just record the max & min each day
SOHChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'','SOH');
SOHds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/soh/%s/%s.%s.D0.SOH.L.%04d.%03d.miniseed','station','network','station','year','jday')
startdate = datenum(2018,3,29);
enddate = datenum(2019,3,26); % last download
%enddate = datenum(2018,4,1)
close all
for stanum = 1:numel(SOHChannelTagList)
    ctag=SOHChannelTagList(stanum);
    channels = {};
    mintime = {};
    maxtime = {};
    miny = {};
    maxy = {};
    for snum=startdate:enddate
        thisdnum = [snum snum+12/24];
        try
            thisw = waveform(SOHds, ctag, snum, snum+1-1/86400);
        catch
            disp('could not read waveform');
            continue;
        end
     
        for channum = 1:numel(thisw) 
            clear ymax ymin indmax indmin tmin tmax
            y = get(thisw(channum), 'data');
            [ymax,indmax] = nanmax(y);
            [ymin,indmin] = nanmin(y);
            try
                dnum = get(thisw(channum),'timevector');
                if numel(y)>=1440 & numel(dnum<1440) % fix for weird sampling rates
                    dnum = snum+(0:1:1439)/1440;
                end
                tmin = dnum(indmin);
                tmax = dnum(indmax);
            catch
%                 size(y)
%                 length(dnum)
%                 indmin
%                 indmax
%                 thisw(c)
%                 close all
%                 plot(y)
                if length(y)>1440
                    [ymax,indmax] = nanmax(y(1:1440));
                    [ymin,indmin] = nanmin(y(1:1440));
                    tmin = dnum(indmin);
                    tmax = dnum(indmax);
                end
            end
            channame = get(thisw(channum),'channel');
            if ~isempty(channels)
                found = false;
                for c=1:numel(channels)
                    if strcmp(channame, channels{c})
                        mintime{c} = [mintime{c} tmin];
                        maxtime{c} = [maxtime{c} tmax];
                        miny{c} = [miny{c} ymin];
                        maxy{c} = [maxy{c} ymax];
                        found = true;
                    end
                end
                if ~found
                    % any channels not found before
                    channels = [channels channame];
                    mintime{channum} = [tmin];
                    maxtime{channum} = [tmax];
                    miny{channum} = [ymin];
                    maxy{channum} = [ymax];
                end
            else
                % for first channel of first file only
                channels = [channels channame];
                mintime{channum} = [tmin];
                maxtime{channum} = [tmax];
                miny{channum} = [ymin];
                maxy{channum} = [ymax];
            end
        end
        save(sprintf('%s_minmaxsoh.mat',ctag.station),'channels', 'mintime', 'maxtime', 'miny',  'maxy');

    end
end

%%
close all
startdate = datenum(2018,3,29);
enddate = datenum(2019,3,26); % last download
for stanum=1:4
    station = sprintf('MTB%d',stanum);

    figure
    load(sprintf('%s_minmaxsoh.mat',station));


    subplot(4,1,1)
    plot(maxtime{5},maxy{5},'b*')
    datetick('x',20)
    xlabel('Date')
    ylabel('# Satellites')
    title(sprintf('Daily max GNSS satellites at %s',station))
    %legend('daily max')
    set(gca,'XLim',[startdate enddate]);

    subplot(4,1,2)
    plot(maxtime{9},maxy{9},'b*')
    % hold on
    % plot(mintime{13},miny{13}/100,'rs')
    datetick('x',20)
    xlabel('Date')
    ylabel('Clock quality (%)')
    title(sprintf('Daily max GNSS clock quality at %s',station))
    %legend('max','min')
    hold off
    set(gca,'YLim',[-1 101],'XLim',[startdate enddate]);

    subplot(4,1,3)
    plot(maxtime{13},maxy{13}/1000,'b*')
    % hold on
    % plot(mintime{13},miny{13}/100,'rs')
    datetick('x',20)
    xlabel('Date')
    ylabel('Voltage (V)')
    title(sprintf('Daily max battery voltage at %s',station))
    %legend('max','min')
    hold off
    set(gca,'YLim',[-0.15 15.15],'XLim',[startdate enddate]);

    subplot(4,1,4)
    ymax = maxy{14};
    ymin = miny{14};
    a = [abs(ymax); abs(ymin)];
    size(a)
    b = max(a,2);
    size(b)
    plot(maxtime{14},b/1000000,'b*')
    % hold on
    % plot(mintime{13},miny{13}/100,'rs')
    datetick('x',20)
    xlabel('Date')
    ylabel('Mass Position (V)')
    title(sprintf('Daily maximum mass position at %s',station))
    %legend('max','min')
    hold off
    set(gca,'YLim',[-0.15 4.15],'XLim',[startdate enddate]);
end

%%
%%
% soh read files into MAT files
clc
SOHChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'','SOH');
%SOHChannelTagList = ChannelTag.array('MV',{'MTB1'},'','SOH');
SOHds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/soh/%s/%s.%s.D0.SOH.L.%04d.%03d.miniseed','station','network','station','year','jday')
startdate = datenum(2018,11,1);
enddate = datenum(2019,3,26); % last download
% startdate = datenum(2018,11,14);
% enddate = datenum(2018,11,16); % last download
close all
for stanum = 1:numel(SOHChannelTagList)
    ctag=SOHChannelTagList(stanum);
    wsoh = [];
    for enum=enddate:-1:startdate
        %try
            thisw = waveform(SOHds, ctag, enum-1, enum);
            thisw = fix_centaur_soh(thisw); % try to fix sample rate, etc.
            thatw = [];
            for c=1:numel(thisw)
                if get(thisw(c),'data_length')>0 % only combine if got data
                    thatw = [thatw thisw(c)];
                end
            end
            if ~isempty(thatw)
                if ~isempty(wsoh)
                    wsoh = combine([wsoh thatw]);
                else 
                    wsoh = thatw;
                end
                wsoh
            else
                disp('no waveform data to combine');
            end
        save(sprintf('soh_from_end_%s_%s_%s.mat',ctag.station,datestr(startdate,'yyyymmdd'),datestr(snum,'yyyymmdd')),'wsoh')
        %end
    end
end
