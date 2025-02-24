import matplotlib.pyplot as plt
import numpy as np
from obspy import Stream


def mulplt2(st, outfile=None, bottomlabel=None, ylabels=None):
    """ Create a plot of a Stream object similar to Seisan's mulplt """
    fh = plt.figure(figsize=(8,12))
    
    # get number of stations
    stations = []
    for tr in st:
        stations.append(tr.stats.station)
    stations = list(set(stations))
    n = len(stations)
    
    # start time as a Unix epoch
    startepoch = st[0].stats.starttime.timestamp
    
    # create empty set of subplot handles - without this any change to one affects all
    axh = []
    
    # loop over all stream objects
    colors = 'kgrb'
    channels = 'ZNEF'
    linewidths = [0.25, 0.1, 0.1, 0.25]
    for i in range(n):
        # add new axes handle for new subplot
        #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
        if i>0:
            #axh.append(plt.subplot(n, 1, i+1, sharex=axh[0]))
            axh.append(fh.add_subplot(n, 1, i+1, sharex=axh[0]))
        else:
            axh.append(fh.add_subplot(n, 1, i+1))
        
        # find all the traces for this station
        this_station = stations[i]
        these_traces = st.copy().select(station=this_station)
        for this_trace in these_traces:
            this_component = this_trace.stats.channel[2]
            line_index = channels.find(this_component)
            #print(this_trace.id, line_index, colors[line_index], linewidths[line_index])
        
            # time vector, t, in seconds since start of record section
            t = np.linspace(this_trace.stats.starttime.timestamp - startepoch,
                this_trace.stats.endtime.timestamp - startepoch,
                this_trace.stats.npts)
            #y = this_trace.data - offset
            y = this_trace.data
            
            # PLOT THE DATA
            print(i, line_index)
            axh[i].plot(t, y, linewidth=linewidths[line_index], color=colors[line_index])
            axh[i].autoscale(enable=True, axis='x', tight=True)
   
        # remove yticks because we will add text showing max and offset values
        #axh[i].yaxis.set_ticks([])

        # remove xticklabels for all but the bottom subplot
        if i < n-1:
            axh[i].xaxis.set_ticklabels([])
        else:
            # for the bottom subplot, also add an xlabel with start time
            if bottomlabel:
                plt.xlabel(bottomlabel)
            else:
                plt.xlabel("Starting at %s" % (st[0].stats.starttime) )

        # default ylabel is station.channel
        if ylabels:
            plt.ylabel(ylabels[i])
        else:
            plt.ylabel(this_station, rotation=90)
            
    # change all font sizes
    plt.rcParams.update({'font.size': 8})
    
    plt.subplots_adjust(wspace=0.1)
    
    # show the figure
    if outfile:
        plt.savefig(outfile, bbox_inches='tight')    
    else:
        plt.show()
    
    
def mulplt(st, bottomlabel='', ylabels=[], MAXPANELS=6):
    """ Create a plot of a Stream object similar to Seisan's mulplt """
    fh = plt.figure()
    n = np.min([MAXPANELS, len(st)])
    
    # start time as a Unix epoch
    startepoch = st[0].stats.starttime.timestamp
    
    # create empty set of subplot handles - without this any change to one affects all
    axh = []
    
    # loop over all stream objects
    for i in range(n):
        # add new axes handle for new subplot
        #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
        axh.append(plt.subplot(n, 1, i+1))
        
        # time vector, t, in seconds since start of record section
        t = np.linspace(st[i].stats.starttime.timestamp - startepoch,
            st[i].stats.endtime.timestamp - startepoch,
            st[i].stats.npts)
            
        # We could detrend, but in case of spikes, subtracting the median may be better
        #st[i].detrend()
        offset = np.median(st[i].data)
        y = st[i].data - offset
        
        # PLOT THE DATA
        axh[i].plot(t, y)
   
        # remove yticks because we will add text showing max and offset values
        axh[i].yaxis.set_ticks([])

        # remove xticklabels for all but the bottom subplot
        if i < n-1:
            axh[i].xaxis.set_ticklabels([])
        else:
            # for the bottom subplot, also add an xlabel with start time
            if bottomlabel=='':
                plt.xlabel("Starting at %s" % (st[0].stats.starttime) )
            else:
                plt.xlabel(bottomlabel)

        # default ylabel is station.channel
        if ylabels==[]:
            plt.ylabel(st[i].stats.station + "." + st[i].stats.channel, rotation=0)
        else:
            plt.ylabel(ylabels[i])

        # explicitly give the maximum amplitude and offset(median)
        plt.text(0, 1, "max=%.1e offset=%.1e" % (np.max(np.abs(y)), offset),
            horizontalalignment='left',
            verticalalignment='top',transform=axh[i].transAxes)
            
    # change all font sizes
    plt.rcParams.update({'font.size': 8})
    
    # show the figure
    plt.show()
    #st.mulplt = types.MethodType(mulplt,st)    
    
    return fh, axh
   
   
    
def plot_stream_types(st, eventdir, maxchannels=10): 
    # assumes input traces are in velocity or pressure units and cleaned
      
    # velocity, displacement, acceleration seismograms
    stV = st.copy().select(channel='[ESBH]H?') 
    if len(stV)>maxchannels:
        stV = stV.select(channel="[ESBH]HZ") # just use vertical components then
        if len(stV)>maxchannels:
            stV=stV[0:maxchannels]
    stD = stV.copy().integrate()   
    stA = stV.copy().differentiate()

    # Infrasound data
    stP = st.copy().select(channel="[ESBH]D?")
    if len(stP)>maxchannels:
        stP=stP[0:maxchannels]
        
    # Plot displacement seismogram
    if stD:
        dpngfile = os.path.join(eventdir, 'seismogram_D.png')
        stD.plot(equal_scale=False, outfile=dpngfile)    
    
    # Plot velocity seismogram
    if stV:
        vpngfile = os.path.join(eventdir, 'seismogram_V.png')
        stV.plot(equal_scale=False, outfile=vpngfile)
         
    # Plot acceleration seismogram
    if stA:
        apngfile = os.path.join(eventdir, 'seismogram_A.png')
        stA.plot(equal_scale=False, outfile=apngfile)
        
    # Plot pressure acoustograms
    if stP:
        ppngfile = os.path.join(eventdir, 'seismogram_P.png')
        stP.plot(equal_scale=False, outfile=ppngfile)    

# from librockets

from obspy.signal.invsim import cosine_taper
from scipy.signal import hilbert, convolve
def get_envelope(tr, seconds=1):
    envelope = np.abs(hilbert(tr.data))
    # Smooth the envelope using a moving average filter
    window_size = int(tr.stats.sampling_rate*seconds)  # Adjust the window size based on your needs
    window = np.ones(window_size) / window_size  # Simple moving average kernel
    envelope = convolve(envelope, window, mode='same')   
    return envelope 

def plot_envelope2(st, window_size=1.0, percentile=99, outfile=None, units=None):
    # Plot the Nth percentile for each 1-second window
    plt.figure(figsize=(10, 6))
    colors = ['red', 'blue', 'green', 'orange', 'black', 'grey', 'purple', 'cyan']
    for i, tr in enumerate(st):

        # Calculate the number of samples per window (1 second)
        samples_per_window = int(window_size * tr.stats.sampling_rate)

        # Get the data and time values
        data = abs(tr.data)
        times = tr.times()

        # Reshape the data into windows (each row is a 1-second window)
        num_windows = len(data) // samples_per_window  # Number of full windows
        reshaped_data = data[:num_windows * samples_per_window].reshape((num_windows, samples_per_window))

        # Calculate the 99th percentile for each window along axis 1 (columns)
        percentiles = np.nanpercentile(reshaped_data, percentile, axis=1)

        # Create the corresponding time values for each window (mid-point of each window)
        window_times = times[:num_windows * samples_per_window:samples_per_window] + window_size / 2


        plt.plot(window_times, percentiles, label=tr.id, color=colors[i % len(colors)], lw=1)
    plt.title(f"{percentile}th Percentile in {window_size}-Second windows")
    plt.xlabel(f"Time (s) from {st[0].stats.starttime}")
    if units:
        plt.ylabel(units)
    plt.grid(True)
    plt.legend()
    if outfile:
        plt.savefig(fname=outfile)
    else:
        plt.show()

def plot_seismograms(st, outfile=None, bottomlabel=None, ylabels=None, units=None, channels='ZNE'):
    """ Create a plot of a Stream object similar to Seisan's mulplt """
    fh = plt.figure(figsize=(8,12))

    from cycler import cycler

    # Define a list of colors you want to cycle through
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # Set the color cycle in matplotlib
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)    
    
    # get number of stations
    stations = []
    for tr in st:
        stations.append(tr.stats.station)
    stations = list(set(stations))
    n = len(stations)
    
    # start time as a Unix epoch
    startepoch = st[0].stats.starttime.timestamp
    
    # create empty set of subplot handles - without this any change to one affects all
    axh = []
    
    # loop over all stream objects
    #colors = ['black', 'blue', 'green']
    
    linewidths = [0.25, 0.1, 0.1, 0.25]
    for i in range(n):
        # add new axes handle for new subplot
        #axh.append(plt.subplot(n, 1, i+1, sharex=ax))
        if i>0:
            #axh.append(plt.subplot(n, 1, i+1, sharex=axh[0]))
            axh.append(fh.add_subplot(n, 1, i+1, sharex=axh[0]))
        else:
            axh.append(fh.add_subplot(n, 1, i+1))
        
        # find all the traces for this station
        this_station = stations[i]
        these_traces = st.copy().select(station=this_station)
        all_ys = []
        lw = 0.5
        if len(these_traces)==1:
            lw = 2
        for this_trace in these_traces:
            this_component = this_trace.stats.channel[2]
            line_index = channels.find(this_component)
            if line_index>-1:
                #print(this_trace.id, line_index, colors[line_index], linewidths[line_index])
            
                # time vector, t, in seconds since start of record section
                t = np.linspace(this_trace.stats.starttime.timestamp - startepoch,
                    this_trace.stats.endtime.timestamp - startepoch,
                    this_trace.stats.npts)
                #y = this_trace.data - offset
                y = get_envelope(this_trace, seconds=0.1)
                all_ys.append(y)
                
                # PLOT THE DATA
                axh[i].plot(t, y, lw=lw, color=colors[line_index], label=this_component)
                axh[i].autoscale(enable=True, axis='x', tight=True)
        if len(these_traces)==3 and len(all_ys)==3:
                vector_amplitude = np.sqrt(all_ys[0]**2 + all_ys[1]**2 + all_ys[2]**2)
                axh[i].plot(t, vector_amplitude, color='red', label='vector', lw=2)
                axh[i].autoscale(enable=True, axis='x', tight=True)  
        elif len(these_traces)>1:
                # Stack the data from each trace
                data = np.array([y for y in all_ys])

                # Compute the median across traces (axis 0 means median across the first dimension - i.e., along traces)
                median_data = np.nanmedian(data, axis=0)
                axh[i].plot(t, median_data, color='red', label='median', lw=2)
                axh[i].autoscale(enable=True, axis='x', tight=True)                         
        plt.grid()
        plt.legend()
        

   
        # remove yticks because we will add text showing max and offset values
        #axh[i].yaxis.set_ticks([])
        '''
        # remove xticklabels for all but the bottom subplot
        if i < n-1:
            axh[i].xaxis.set_ticklabels([])
        else:
            # for the bottom subplot, also add an xlabel with start time
            if bottomlabel:
                plt.xlabel(bottomlabel)
            else:
                plt.xlabel("Starting at %s" % (st[0].stats.starttime) )
        '''

        # default ylabel is station.channel
        ylabelstr = this_station + '\n' + units
        if ylabels:
            ylabelstr = ylabels[i]
        plt.ylabel(ylabelstr, rotation=90)

        plt.xlabel(f'Seconds from {st[0].stats.starttime}')

    # change all font sizes
    plt.rcParams.update({'font.size': 10})
    
    plt.subplots_adjust(wspace=0.1)

    #plt.suptitle(f'Amplitude from {st[0].stats.starttime}')
    
    # show the figure
    if outfile:
        plt.savefig(outfile, bbox_inches='tight')    
    else:
        plt.show()
    

def floor_minute(timestamp):
    return timestamp.replace(second=0, microsecond=0) 

def ceil_minute(timestamp):
    return (timestamp+60).replace(second=0, microsecond=0) 

def super_stream_plot(st, dpi=100, Nmax=6, rank='ZFODNE123456789', equal_scale=True, outfile=None, figsize=(8,8)):
    fig = plt.figure(figsize=figsize)  # Set a custom figure size (width, height)
    if len(st)>Nmax:
        st2 = obspy.Stream()
        rankpos = 0
        while len(st2)<Nmax and rankpos<len(rank):
            for tr in st.select(channel=f'*{rank[rankpos]}'):
                if len(st2)<Nmax:
                    st2.append(tr)
            rankpos += 1
    else:
        st2 = st.copy()

    height = 250
    if len(st2)>4:
        height=1000/len(st2)
        #height=100
    #size_tuple = (1200, height) * len(st2)
    #size_tuple = (height, 1200) * len(st2)
    size_tuple = (1000, height) * len(st2)
    st2.plot(fig=fig, dpi=dpi, color='b', equal_scale=equal_scale, outfile=outfile);


def downsample_and_plot_1_trace_per_location(stream, target_sampling_rate=1.0):
        # Downsample to 1 Hz if necessary
        print(f'Downsampling traces', end=' ')
        for trace in stream:
            print(trace.id, end=' ')
            if trace.stats.sampling_rate > target_sampling_rate:
                decimation_factor = int(trace.stats.sampling_rate / target_sampling_rate)
                if decimation_factor > 1:
                    trace.decimate(decimation_factor, no_filter=True)
        print('\n')

        # Group traces by station-location and select the ones with the most non-zero samples
        station_location_data = {}
        print('Subsetting channels to plot')
        for trace in stream:
            key = (trace.stats.station, trace.stats.location)
            channel = trace.stats.channel
            nonzero_count = np.count_nonzero(trace.data)

            if key not in station_location_data:
                station_location_data[key] = {"seismic": None, "infrasound": None}

            # Seismic selection: Highest non-zero sample count among ?HZ, ?HN, ?HE
            if channel[1]=='H':
                if (
                    station_location_data[key]["seismic"] is None or
                    nonzero_count > count_nonzero_samples(station_location_data[key]["seismic"])
                ):
                    station_location_data[key]["seismic"] = trace  # Pick the trace with the most non-zero samples

            # Infrasound selection: Highest non-zero sample count among ?DF, ?DG, ?D
            if channel[1]=='D':
                if (
                    station_location_data[key]["infrasound"] is None or
                    nonzero_count > count_nonzero_samples(station_location_data[key]["infrasound"])
                ):
                    station_location_data[key]["infrasound"] = trace  # Pick the trace with the most non-zero samples

        # Create a single stream for all selected traces
        selected_stream = Stream()
        for channels in station_location_data.values():
            if channels["seismic"]:
                selected_stream += channels["seismic"]
            if channels["infrasound"]:
                selected_stream += channels["infrasound"]

        # Plot all selected traces in one figure
        print('Plotting')
        if len(selected_stream) > 0:
            print(f"Plotting {len(selected_stream)} traces in one figure.")
            selected_stream.plot(equal_scale=False)



def plot_station_amplitude_map(st, station0hypfile=None, outfile=None):
    import cartopy.crs as crs
    import cartopy.feature as cf
    from seisan_classes import parse_STATION0HYP, add_station_locations
    if not station0hypfile:
        return
    
    station_locationsDF = parse_STATION0HYP(station0hypfile)
    add_station_locations(st, station_locationsDF)  
    
    # draw Montserrat coastline
    extent = [-62.27, -62.12, 16.67, 16.82]
    #central_lon = np.mean(extent[:2])
    #central_lat = np.mean(extent[2:])
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1,1,1, projection=crs.PlateCarree())
    ax.set_extent(extent, crs=crs.PlateCarree())
    ax.add_feature(cf.COASTLINE)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    
    l = len(st)
    amps = np.zeros((l, 1))
    lons = np.zeros((l, 1))
    lats = np.zeros((l, 1))
    #calibs = []
    textlabels = []
    for i, tr in enumerate(st):
        lons[i] = tr.stats.lon
        lats[i] = tr.stats.lat
        amps[i] = tr.stats.metrics.peakamp
        #calibs.append(tr.stats.calib)
        textlabels.append(tr.stats.station)
    sizes = amps/max(amps) 
    #print('calibs=',calibs)
    g = ax.scatter(x=lons, y=lats,color="red",s=sizes*300,alpha=0.8,transform=crs.PlateCarree());
    g.set_facecolor('none');
    g.set_edgecolor('red');    
    for i,thislabel in enumerate(textlabels):
        ax.text(lons[i],lats[i],thislabel,transform=crs.PlateCarree());
            
    # add the volcano
    ax.scatter(-62.1833326, 16.7166638,  300, marker='*', transform=crs.PlateCarree());
         
    # save or show
    if outfile:
        try:
            plt.savefig(outfile,  bbox_inches='tight');
        except:
            plt.show();
    else:
        plt.show();