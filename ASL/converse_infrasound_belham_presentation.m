type('readme.m');
settings;
eventnum = selectevent(snum, enum);

% Analyze event
analyze_flood(TOPDIR, station, ds, snum(eventnum), enum(eventnum))

%%
function w=analyze_flood(TOPDIR, station, ds, snum, enum)
    % load Belham network data from original Miniseed files, or from MAT files if they exist already
    disp(sprintf('Loading data for event starting at %s',datestr(snum,30)))
    wv = [];
    ChannelTagList = [station.ctags];
    
    if ~exist('calibObjects','var')
        calibObjects = iceweb.usf_calibrations(ChannelTagList); % retrieve calibrations for Trillium and InfraBSU sensors
    end    
    
    for c=1:numel(ChannelTagList)
        ctag = ChannelTagList(c);
        fprintf('Processing %s',ctag.string())
        MATfilename = sprintf('%s/%s.%s.mat',TOPDIR,ctag.string(),datestr(snum(1),30));
        if exist(MATfilename,'file')
            fprintf('- Loading %s\n',MATfilename)
            load(MATfilename,'w');
        else
            fprintf('- Loading from Miniseed\n')
            w = waveform(ds, ChannelTagList(c), snum, enum);
            % Clean data and convert to correct units
            w = clean(w, 0.1); % applies a high-pass 10-second (0.1 Hz) filter
            w = w * calibObjects(c).calib;
            w = set(w,'units',calibObjects(c).units);        
            fprintf('- Saving to %s\n',MATfilename)
            save(MATfilename,'w') % create MAT files for faster loading next time
        end
        wv = [wv w];
        clear w
    end
    disp('- done');
    
    wv_orig = wv;
    distances = [station.distance];
    choice = -1;
    while choice < 13
        choice=menu('Menu', 'Extract pulse', 'Select infrasound only', 'Select vertical seismic only', 'Pick onsets by hand', 'Spectra', 'Spectrograms', 'Plot RSAM', 'Record section', 'Panel plot', 'Back to original traces', 'Slope plot', 'Plot in frequency bands','quit');
        close all
        switch choice
            case 1, [x y]=extractpulse(wv); dnum = x/86400 + get(wv(1),'start'); wv = extract(wv, 'time', min(dnum), max(dnum));
            case 2, wv=selectchannels(wv,'HDF');
            case 3, wv=selectchannels(wv,'HHZ');
            case 4, [wv,detObj]=pickonsets(wv,'manual', distances);
            %case 5, [wv,detObj]=pickonsets(wv,'automatic', distances);
            case 5, plotspectra(wv);
            case 6, plotspectrograms(wv);
            case 7, plotrsam(wv);
            case 8, if numel(distances) ~= numel(wv); distances = sort(repmat(distances,1,4)); end; recordsection(w, distances);
            case 9, plot_panels(wv);
            case 10, wv = wv_orig;
            case 11, slopeplot(distances, elevations, stations);
            case 12, plotbands(wv);
            case 13, break;
        end
    end
end













% %% Sound file
% waveform2sound(wv_clean(4), 120, 'laharsound1.wav')








