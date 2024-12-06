% extract a few hours around earthquake 
% https://earthquake.usgs.gov/earthquakes/eventpage/us1000jkrq/executive
% M 6.1 - 7km NW of El Dovio, Colombia
% 2019-03-23 19:21:16 (UTC)4.560°N 76.280°W113.3 km depth
ds = datasource('miniseed', '~/Dropbox/MY_BRIEFCASE/montserrat_projects/Belham_Valley_Project/belhamdownload3/%04d/%s/%s/%s.D/%s.%s.%s.%s.D.%04d.%03d.miniseed','year','network','station','channel','network','station','location','channel','year','jday')
ChannelTagList = ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHE');
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHN')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'00','HHZ')];
ChannelTagList = [ChannelTagList ChannelTag.array('MV',{'MTB1';'MTB2';'MTB3';'MTB4'},'10','HDF')];
snum = datenum(2019,3,23,19,0,0);
enum = datenum(2019,3,23,20,0,0); 
w = waveform(ds, ChannelTagList, snum, enum);
calibObjects = iceweb.usf_calibrations(ChannelTagList);
w2 = iceweb.apply_calib(w, calibObjects);
save eq.mat
%plot_panels(w([1:3 5:7 9:11 13:15])) % ignore infrasound
%%
stations = unique({ChannelTagList.station});
for stanum = 1:numel(stations)
    figure
    firstchan = (stanum-1)*4;
    ax=[];
    for c=1:3
        ax = subplot(3,1,c);
        plot(w(firstchan+c),'axeshandle',ax);
    end
end
%%
s1 = get(w(1),'data');
s2 = get(w(2),'data');
%
s3 = s2(2e5:6e5);
[acor,lag] = xcorr(s3,s1);
%
[~,I] = max(abs(acor));
timeDiff = lag(I);
figure
subplot(311); plot(s1); title('s1');
subplot(312); plot(s3); title('s3');
subplot(313); plot(lag,acor);
%
timeDiff
%%
close all

% set the STA/LTA detector
sta_seconds = 5; % STA time window 0.7 seconds
lta_seconds = 50; % LTA time window 7 seconds
thresh_on = 3; % Event triggers "ON" with STA/LTA ratio exceeds 3
thresh_off = 1.5; % Event triggers "OFF" when STA/LTA ratio drops below 1.5
minimum_event_duration_seconds = 20.0; % Trigger must be on at least 2 secs
pre_trigger_seconds = 0; % Do not pad before trigger
post_trigger_seconds = 0; % Do not pad after trigger
event_detection_params = [sta_seconds lta_seconds thresh_on thresh_off ...
    minimum_event_duration_seconds];

%%
% run the STA/LTA detector to fix the clock drift on MTB1, MTB2, MTB4
% relative so far - I haven't kept the time of MTB3 fixed, but I should
% because it is the only station with a good GPS clock
close all
detectiontime = NaN(numel(w),1);
for wavnum=1:numel(w)
    [detObj,sta,lta,sta_to_lta] = Detection.sta_lta(w(wavnum), 'edp', event_detection_params, ...
        'lta_mode', 'frozen');
    if strcmp(class(detObj),'Detection')
        detectiontime(wavnum) = detObj.time(1);
    end
end
%%
mtb1_detectiontime = nanmedian(detectiontime(1:4));
mtb2_detectiontime = nanmedian(detectiontime(5:8));
mtb3_detectiontime = nanmedian(detectiontime(9:12));
mtb4_detectiontime = nanmedian(detectiontime(13:16));
mtb1_shifttime = mtb3_detectiontime - mtb1_detectiontime;
mtb2_shifttime = mtb3_detectiontime - mtb2_detectiontime;
mtb3_shifttime = mtb3_detectiontime - mtb3_detectiontime;
mtb4_shifttime = mtb3_detectiontime - mtb4_detectiontime;
w2 = [];
for c=1:4
    oldstarttime = get(w(c),'start');
    newstarttime = oldstarttime + mtb1_shifttime;
    w2 = [w2 set(w(c),'start', newstarttime) ];
end
for c=5:8
    oldstarttime = get(w(c),'start');
    newstarttime = oldstarttime + mtb2_shifttime;
    w2 = [w2 set(w(c),'start', newstarttime) ];
end
w2(9:12) = w(9:12);
for c=13:16
    oldstarttime = get(w(c),'start');
    newstarttime = oldstarttime + mtb3_shifttime;
    w2 = [w2 set(w(c),'start', newstarttime) ];
end
plot_panels(w2);

    