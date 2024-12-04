import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import struct
import fnmatch
import obspy

class RSAM:

    def __init__(self, dataframes=None, stream=None, sampling_interval=60.0, filter=[0.5, 18.0], bands = {'VLP': [0.02, 0.2], 'LP':[0.5, 4.0], 'VT':[4.0, 18.0]}):
        ''' Create an RSAM object 
        
            Optional name-value pairs:
                dataframes: Creates an RSAM object using these dataframes. Used by downsample() method, for example. Default: None.
                stream: Creates as RSAM object from this ObsPy.Stream object.
                sampling_interval: Compute RSAM data using this sampling interval (in seconds). Default: 60
                filter: list of two floats, representing fmin and fmax. Default: [0.5, 18.0]. Set to None if no filter wanted.
                bands: a dictionary of filter bands and corresponding column names. Default: {'VLP': [0.02, 0.2], 'LP':[0.5, 4.0], 
                    'VT':[4.0, 18.0]}. For example, the default setting creates 3 additional columns for each DataFrame called 
                    'VLP', 'LP', and 'VT', which contain the mean value for each sampling_interval within the specified filter band
                    (e.g. 0.02-0.2 Hz for VLP). If 'LP' and 'VT' are in this dictionary, an extra column called 'fratio' will also 
                    be computed, which is the log2 of the ratio of the 'VT' column to the 'LP' column, following the definition of
                    frequency ratio by Rodgers et al. (2015).
        '''
        self.dataframes = {} 
        #self.trace_ids = []

        if isinstance(dataframes, dict):
            good_dataframes = {}
            for id, df in dataframes.items():
                if isinstance(df, pd.DataFrame):
                    good_dataframes[id]=df
            if len(good_dataframes)>0:
                self.dataframes = good_dataframes
                print('dataframes found. ignoring other arguments.')
                return
            else:
                print('no valid dataframes found')

        if not isinstance(stream, obspy.core.Stream):
            # empty RSAM object
            print('creating blank RSAM object')
            return

        if len(stream)>0:
            if stream[0].stats.sampling_rate == 1/sampling_interval:
                # no downsampling to do
                for tr in stream:
                    df = pd.DataFrame()
                    df['time'] = pd.Series(tr.times('timestamp'))
                    df['mean'] = pd.Series(tr.data) 
                    self.dataframes[tr.id] = df
                    #self.trace_ids.append(tr.id)  
                return 
            elif stream[0].stats.sampling_rate < 1/sampling_interval:
                print('error: cannot compute RSAM for a Stream with a tr.stats.delta bigger than requested sampling interval')
                return
            
        for tr in stream:
            if tr.stats.npts < tr.stats.sampling_rate * sampling_interval:
                print('Not enough samples for ',tr.id,'. Skipping.')
                continue
            #print(tr.id, 'absolute=',absolute)
            df = pd.DataFrame()
            
            t = tr.times('timestamp') # Unix epoch time
            sampling_rate = tr.stats.sampling_rate
            t = self.__reshape_trace_data(t, sampling_rate, sampling_interval)
            df['time'] = pd.Series(np.nanmin(t,axis=1))

            if filter:
                if tr.stats.sampling_rate<filter[1]*2.5:
                    #print(f"{tr}: bad sampling rate. Skipping.")
                    continue
                tr2 = tr.copy()
                tr2.detrend('demean')
                tr2.filter('bandpass', freqmin=filter[0], freqmax=filter[1], corners=2)
                y = self.__reshape_trace_data(np.absolute(tr2.data), sampling_rate, sampling_interval)
            else:
                y = self.__reshape_trace_data(np.absolute(tr.data), sampling_rate, sampling_interval)

            df['min'] = pd.Series(np.nanmin(y,axis=1))   
            df['mean'] = pd.Series(np.nanmean(y,axis=1)) 
            df['max'] = pd.Series(np.nanmax(y,axis=1))
            df['median'] = pd.Series(np.nanmedian(y,axis=1))
            df['std'] = pd.Series(np.nanstd(y,axis=1))

            if bands:
                for key in bands:
                    tr2 = tr.copy()
                    [flow, fhigh] = bands[key]
                    tr2.filter('bandpass', freqmin=flow, freqmax=fhigh, corners=2)
                    y = self.__reshape_trace_data(abs(tr2.data), sampling_rate, sampling_interval)
                    df[key] = pd.Series(np.nanmean(y,axis=1))
                if 'LP' in bands and 'VT' in bands:
                    df['fratio'] = np.log2(df['VT']/df['LP'])
                
            self.dataframes[tr.id] = df
            #self.trace_ids.append(tr.id)
    
    
    def copy(self):
        ''' make a full copy of an RSAM object and return it '''
        selfcopy = RSAM(stream=obspy.core.Stream())
        selfcopy.trace_ids = self.trace_ids
        selfcopy.dataframes = self.dataframes.copy()
        return selfcopy

    def downsample(self, new_sampling_interval=3600):
        ''' downsample an RSAM object to a larger sampling interval(e.g. from 1 minute to 1 hour). Returns a new RSAM object.
         
            Optional name-value pair:
                new_sampling_interval: the new sampling interval (in seconds) to downsample to. Default: 3600
        '''

        dataframes = {}
        for id in self.dataframes:
            df = self.dataframes[id]
            df['date'] = pd.to_datetime(df['time'], unit='s')
            old_sampling_interval = self.__get_sampling_interval(df)
            if new_sampling_interval > old_sampling_interval:
                freq = '%.0fmin' % (new_sampling_interval/60)
                new_df = df.groupby(pd.Grouper(key='date', freq=freq)).mean()
                new_df.reset_index(drop=True)
                dataframes[id] = new_df
            else:
                print('Cannot downsample to a smaller sampling interval')
        return self.__class__(dataframes=dataframes) 
        
    def drop(self, id):
        if id in self.__get_trace_ids():
            del self.dataframes[id]

    def plot(self, metrics=['mean'], kind='stream', logy=False, equal_scale=False):
        ''' plot an RSAM object 

            Optional name-value pairs:
                metrics: The columns of each RSAM DataFrame to plot. Can be one (scalar), or many (a list)
                         If metrics='bands', this is shorthand for metrics=['VLP', 'LP', 'VT', 'specratio']
                         Default: metrics='mean'
                kind:    The kind of plot to make. kind='stream' (default) will convert each of the request 
                         DataFrame columns into an ObsPy.Stream object, and then use the ObsPy.Stream.plot() method.
                         kind='line' will render plots directly using matplotlib.pyplot, with all metrics requested 
                         on a single plot.
                logy:    In combination with kind='line', will make the y-axis logarithmic. No effect if kind='stream'.
                equal_scale: If True, y-axes for each plot will have same limits. Default: False.
        
        '''
        self.__remove_empty()
        if isinstance(metrics, str):
            metrics = [metrics]
        if kind == 'stream':
            if metrics == ['bands']:
                metrics = ['VLP', 'LP', 'VT', 'specratio']
            for m in metrics:
                print('METRIC: ',m)
                st = self.to_stream(metric=m)
                st.plot(equal_scale=equal_scale);
            return
        for key in self.dataframes:
            df = self.dataframes[key]
            this_df = df.copy()
            this_df['time'] = pd.to_datetime(df['time'], unit='s')
            if metrics == ['bands']:
                # plot f-bands only
                if not 'VLP' in this_df.columns:
                    print('no frequency bands data for ',key)
                    continue
                ph2 = this_df.plot(x='time', y=['VLP', 'LP', 'VT'], kind='line', title=f"{key}, f-bands", logy=log, rot=45)
                plt.show()
            else:
                for m in metrics:
                    got_all_metrics = True
                    if not m in this_df.columns:
                        print(f'no {m} column for {key}')
                        got_all_metrics = False
                if not got_all_metrics:
                    continue
                if kind == 'line':
                    ph = this_df.plot(x='time', y=metrics, kind=kind, title=key, logy=logy, rot=45)
                elif kind  == 'scatter':
                    fh, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
                    for i, m in enumerate(metrics):
                        this_df.plot(x='time', y=m, kind=kind, ax=ax[i], title=key, logy=logy, rot=45)
                plt.show()
            plt.close('all')

    @classmethod
    def readRSAMfile(classref, startt, endt, RSAM_TOP, trace_ids=None, sampling_interval=60, ext='pickle'):
        ''' read one or many RSAM files from folder specified by RSAM_TOP for date/time range specified by startt, endt
            return corresponding RSAM object

            startt and endt must be ObsPy.UTCDateTime data types

            Optional name-value pairs:
                trace_ids (list): only load RSAM files corresponding to these trace IDs.
                sampling_interval (int): seconds of raw seismic data corresponding to each RSAM sample. Default: 60
                ext (str): should be 'csv' or 'pickle' (default). Indicates what type of file format to open.

        '''
        #self = classref() # blank RSAM object
        dataframes = {}

        if not trace_ids: # make a list of possible trace_ids, regardless of year
            rsamfiles = glob.glob(os.path.join(RSAM_TOP,'RSAM_*_[0-9][0-9][0-9][0-9]_%ds.%s' % (sampling_interval, ext )))
            trace_ids = []
            for rsamfile in rsamfiles:
                parts = rsamfile.split('_')
                trace_ids.append(parts[-3])
        
        for id in trace_ids:
            df_list = []
            for yyyy in range(startt.year, endt.year+1):
                rsamfile = os.path.join(RSAM_TOP,'RSAM_%s_%4d_%ds.%s' % (id, yyyy, sampling_interval, ext))
                if os.path.isfile(rsamfile):
                    print('Reading ',rsamfile)
                    if ext=='csv':
                        df = pd.read_csv(rsamfile, index_col=False)
                    elif ext=='pickle':
                        df = pd.read_pickle(rsamfile)
                    if df.empty:
                        continue
                    df['pddatetime'] = pd.to_datetime(df['time'], unit='s')
                    # construct Boolean mask
                    mask = df['pddatetime'].between(startt.isoformat(), endt.isoformat())
                    # apply Boolean mask
                    subset_df = df[mask]
                    subset_df = subset_df.drop(columns=['pddatetime'])
                    df_list.append(subset_df)
            if len(df_list)==1:
                dataframes[id] = df_list[0]
            elif len(df_list)>1:
                dataframes[id] = pd.concat(df_list)
            #self.trace_ids.append(id)   
                
        rsamObj = classref(dataframes=dataframes) # create RSAM object         
        return rsamObj

    @classmethod
    def readRSAMbinary(classref, RSAM_DIR, station, stime, etime):
        ''' read one (or many if station is a list) RSAM binary file(s) recorded by the original RSAM system
            return corresponding RSAM object '''
        st = obspy.core.Stream()

        if isinstance(station, list):
            for this_station in station:
                tr = classref.readRSAMbinary(RSAM_DIR, this_station, stime, etime)
                if tr.data.size - np.count_nonzero(np.isnan(tr.data)): # throw away Trace objects with only NaNs
                    st.append(tr)
            rsamObj = classref(stream=st, sampling_interval = 1/st[0].stats.sampling_rate)
            return rsamObj
        else:

            for year in range(stime.year, etime.year+1):
            
                daysPerYear = 365
                if year % 4 == 0:
                    daysPerYear += 1
                    
                RSAMbinaryFile = os.path.join(RSAM_DIR, f"{station}{year}.DAT")
                
                values = []
                if os.path.isfile(RSAMbinaryFile):
                    print('Reading ',RSAMbinaryFile)
        
                    # read the whole file
                    f = open(RSAMbinaryFile, mode="rb")
                    f.seek(4*1440) # 1 day header
                    for day in range(daysPerYear):
                        for minute in range(60 * 24):
                            v = struct.unpack('f', f.read(4))[0]
                            values.append(v)
                            #print(type(v)) 
                            #print(v)
                    f.close()
        
                    # convert to Trace object
                    tr = obspy.Trace(data=np.array(values))
                    tr.stats.starttime = obspy.core.UTCDateTime(year, 1, 1, 0, 0, 0)
                    tr.id=f'MV.{station}..EHZ'
                    tr.stats.sampling_rate=1/60
                    tr.data[tr.data == -998.0] = np.nan
        
                    # Trim based on stime & etime & append to Stream
                    tr.trim(starttime=stime, endtime=etime)
                    st.append(tr)
                else:
                    print(f"{RSAMbinaryFile} not found")
        
            return st.merge(method=0, fill_value=np.nan)[0]


    def select(self, network=None, station=None, location=None, channel=None,
               sampling_interval=None, npts=None, component=None, id=None,
               inventory=None):
        """
        Return new RSAM object only with DataFrames that match the given
        criteria (e.g. all DataFrames with ``channel="BHZ"``).

        Alternatively, DataFrames can be selected based on the content of an
        :class:`~obspy.core.inventory.inventory.Inventory` object: DataFrame will
        be selected if the inventory contains a matching channel active at the
        DataFrame start time.

        based on obspy.Stream.select()

        .. rubric:: Examples

        >>> rsamObj2 = rsamObj.select(station="R*")
        >>> rsamObj2 = rsamObj.select(id="BW.RJOB..EHZ")
        >>> rsamObj2 = rsamObj.select(component="Z")
        >>> rsamObj2 = rsamObj.select(network="CZ")
        >>> rsamObj2 = rsamObj.select(inventory=inv)
    
        All keyword arguments except for ``component`` are tested directly
        against the respective entry in the :class:`~obspy.core.trace.Stats`
        dictionary.

        If a string for ``component`` is given (should be a single letter) it
        is tested against the last letter of the ``Trace.stats.channel`` entry.

        Alternatively, ``channel`` may have the last one or two letters
        wildcarded (e.g. ``channel="EH*"``) to select all components with a
        common band/instrument code.

        All other selection criteria that accept strings (network, station,
        location) may also contain Unix style wildcards (``*``, ``?``, ...).
        """
        if inventory is None:
            dataframes = self.dataframes
        else:
            trace_ids = []
            start_dates = []
            end_dates = []
            for net in inventory.networks:
                for sta in net.stations:
                    for chan in sta.channels:
                        id = '.'.join((net.code, sta.code,
                                       chan.location_code, chan.code))
                        trace_ids.append(id)
                        start_dates.append(chan.start_date)
                        end_dates.append(chan.end_date)
            dataframes = {}
            for thisid, thisdf in self.dataframes.items():
                idx = 0
                while True:
                    try:
                        idx = trace_ids.index(thisid, idx)
                        start_date = start_dates[idx]
                        end_date = end_dates[idx]
                        idx += 1
                        if start_date is not None and\
                                self.__get_starttime(thisdf) < start_date:
                            continue
                        if end_date is not None and\
                                self.__get_endtime(thisdf) > end_date:
                            continue
                        dataframes[thisid]=thisdf
                    except ValueError:
                        break
        dataframes_after_inventory_filter = dataframes

        # make given component letter uppercase (if e.g. "z" is given)
        if component is not None and channel is not None:
            component = component.upper()
            channel = channel.upper()
            if (channel[-1:] not in "?*" and component not in "?*" and
                    component != channel[-1:]):
                msg = "Selection criteria for channel and component are " + \
                      "mutually exclusive!"
                raise ValueError(msg)

        # For st.select(id=) without wildcards, use a quicker comparison mode:
        quick_check = False
        quick_check_possible = (id is not None
                                and sampling_rate is None and npts is None
                                and network is None and station is None
                                and location is None and channel is None
                                and component is None)
        if quick_check_possible:
            no_wildcards = not any(['?' in id or '*' in id or '[' in id])
            if no_wildcards:
                quick_check = True
                [net, sta, loc, chan] = id.upper().split('.')

        dataframes = {}
        for thisid, thisdf in dataframes_after_inventory_filter.items():
            [thisnet, thissta, thisloc, thischan] = thisid.upper().split('.')
            if quick_check:
                if (thisnet.upper() == net
                        and thissta.upper() == sta
                        and thisloc.upper() == loc
                        and thischan.upper() == chan):
                    dataframes.append(thisdf)
                continue
            # skip trace if any given criterion is not matched
            if id and not fnmatch.fnmatch(thisid.upper(), id.upper()):
                continue
            if network is not None:
                if not fnmatch.fnmatch(thisnet.upper(),
                                       network.upper()):
                    continue
            if station is not None:
                if not fnmatch.fnmatch(thissta.upper(),
                                       station.upper()):
                    continue
            if location is not None:
                if not fnmatch.fnmatch(thisloc.upper(),
                                       location.upper()):
                    continue
            if channel is not None:
                if not fnmatch.fnmatch(thischan.upper(),
                                       channel.upper()):
                    continue
            if sampling_interval is not None:
                if float(sampling_interval) != self.__get_sampling_interval(thisdf):
                    continue
            if npts is not None and int(npts) != self.__get_npts(thisdf):
                continue
            if component is not None:
                if not fnmatch.fnmatch(thischan[-1].upper(),
                                       component.upper()):
                    continue
            dataframes[thisid]=thisdf
        return self.__class__(dataframes=dataframes)       
     
    def to_stream(self, metric='mean'):
        ''' Convert one column (specified by metric) of each DataFrame in an RSAM object to an obspy.Trace, 
            return an ObsPy.Stream that is the combination of all of these Trace objects'''
        st = obspy.core.Stream()
        for key in self.dataframes:
            #print(key)
            df = self.dataframes[key]
            if metric in df.columns:
                dataSeries = df[metric]
                tr = obspy.core.Trace(data=np.array(dataSeries))
                #timeSeries = pd.to_datetime(df['time'], unit='s')
                tr.stats.delta = self.__get_sampling_interval(df)
                tr.stats.starttime = obspy.core.UTCDateTime(df.iloc[0]['time'])
                tr.id = key
                if tr.data.size - np.count_nonzero(np.isnan(tr.data)):
                    st.append(tr)
                
        return st
    
    def trim(self, starttime=None, endtime=None, pad=False, keep_empty=False, fill_value=None):
        ''' trim RSAM object based on starttime and endtime. Both must be of type obspy.UTCDateTime 

            based on obspy.Stream.trim()

            keep_empty=True will retain dataframes that are either blank, or full of NaN's or 0's.
                       Default: False
            
            Note:
            - pad option and fill_value option not yet implemented
        '''
        if pad:
            print('pad option not yet supported')
            return
        if fill_value:
            print('fill_value option not yet supported')
            return
        if not starttime or not endtime:
            print('starttime and endtime required as ObsPy.UTCDateTime')
            return
        for id in self.dataframes:
            df = self.dataframes[id]
            mask = (df['time']  >= starttime.timestamp ) & (df['time'] <= endtime.timestamp )
            self.dataframes[id] = df.loc[mask]
        if not keep_empty:
            self.__remove_empty()
        
    def write(self, RSAM_TOP, ext='pickle'):
        ''' Write RSAM object to CSV or Pickle files (one per net.sta.loc.chan, per year) into folder specified by RSAM_TOP

            Optional name-value pairs:
                ext: Should be 'csv' or 'pickle' (default). Specifies what file format to save each DataFrame.   
        '''
        if not os.path.isdir(RSAM_TOP):
            os.makedirs(RSAM_TOP)
        for key in self.dataframes:
            df = self.dataframes[key]
            if df.empty:
                continue
            starttime = df.iloc[0]['time']
            yyyy = obspy.core.UTCDateTime(starttime).year
            rsamfile = os.path.join(RSAM_TOP,'RSAM_%s_%4d_%ds.%s' % (key, yyyy, self.__get_sampling_interval(df), ext))

            #print(f"Saving to {rsamfile}")
            if os.path.isfile(rsamfile):
                if ext=='csv':
                    original_df = pd.read_csv(rsamfile)
                elif ext=='pickle':
                    original_df = pd.read_pickle(rsamfile)
                # SCAFFOLD: should check if RSAM data already exist in file for the DataFrame time range.
                # Currently just combining and dropping duplicates without much thought. maybe ok?
                combined_df = pd.concat([original_df, df], ignore_index=True)
                combined_df = combined_df.drop_duplicates(subset=['time'], keep='last') # overwrite duplicate data
                print(f'Modifying {rsamfile}')
                if ext=='csv':
                    combined_df.to_csv(rsamfile, index=False)
                elif ext=='pickle':
                    combined_df.to_pickle(rsamfile)
            else:
                # SCAFFOLD: do i need to create a blank file here for whole year? probably not because smooth() is date-aware
                print(f'Writing {rsamfile}')
                if ext=='csv':
                    df.to_csv(rsamfile, index=False)
                elif ext=='pickle':
                    df.to_pickle(rsamfile)


    @staticmethod
    def __get_endtime(df):
        ''' return the end time of an RSAM dataframe as an ObsPy UTCDateTime'''
        return obspy.core.UTCDateTime(df.iloc[-1]['time'])
    

    @staticmethod
    def __get_npts(df):
        ''' return the number of rows of an RSAM dataframe'''
        return len(df)

    @staticmethod
    def __get_sampling_interval(df):
        ''' return the sampling interval of an RSAM dataframe in seconds '''
        return df.iloc[1]['time'] - df.iloc[0]['time']       

    @staticmethod
    def __get_starttime(df):
        ''' return the start time of an RSAM dataframe as an ObsPy UTCDateTime'''
        return obspy.core.UTCDateTime(df.iloc[0]['time'])
    
    def __get_trace_ids(self):
        return [id for id in self.dataframes]

            
    def __remove_empty(self):
        ''' remove empty dataframes from an RSAM object - these are net.sta.loc.chan for which there are no non-zero data '''
        #print('removing empty dataframes')
        dfs_dict = self.dataframes.copy() # store this so we can delete during loop, otherwise complains about deleting during iteration
        for id in self.dataframes:
            #print(id, self.dataframes[id]['mean'])
            if len(self.dataframes[id])==0 or (self.dataframes[id]['mean'] == 0).all() or (pd.isna(self.dataframes[id]['mean'])).all():
                del dfs_dict[id]
                self.trace_ids.remove(id)
        self.dataframes = dfs_dict

    @staticmethod
    def __reshape_trace_data(x, sampling_rate, sampling_interval):
        ''' reshape data vector from 1-D to 2-D to support vectorized loop for RSAM computation '''
        # reshape the data vector into an array, so we can take advantage of np.mean()
        x = np.absolute(x)
        s = np.size(x) # find the size of the data vector
        nc = int(sampling_rate * sampling_interval) # number of columns
        nr = int(s / nc) # number of rows
        x = x[0:nr*nc] # cut off any trailing samples
        y = x.reshape((nr, nc))
        return y

    def __str__(self):
        contents=""
        for i, trid in enumerate(self.dataframes):
            df = self.dataframes[trid]
            if i==0:
                contents += f"Metrics: {','.join(df.columns[1:])}" + '\n'
                contents += f"Sampling Interval={self.__get_sampling_interval(df)} s" + '\n\n'
            startt = self.__get_starttime(df)
            endt = self.__get_endtime(df)        
            contents += f"{trid}: {startt.isoformat()} to {endt.isoformat()}"
            contents += "\n"
        return contents   
    
class Dr(RSAM): # inherit from RSAM module
    # how will this be different?
    # - accept an inventory
    # - remove instrument response and integrate to displacement
    # - compute distances
    # - correct for geometrical spreading (body or surface waves)
    # - optionally correct for inelastic attenuation
    # - make an IceWeb-style Dr plot
    # 
    # Design question: - just correct for geometrical spreading and attenuation when plotting?
    #                  - this would mean saving Displacement, but not Reduced Displacement
    #                  - could write an intermediate class
    #
    # What methods would change, and stay same:
    # 

class ReducedDisplacementObj(): # SCAFFOLD: Modify this so it inherits from RSAM and modifies just a few methods, and adds some
    # Difference from RSAM is 
    # (1) Correct to displacement first 
    # (2) Filter, 
    # (3) Apply geometric spreading 
    # (4) Correct for Q
    def __init__(self, st=None, inv=None, sampling_interval=60.0, verbose=False, metric='median', \
            freqmin=0.5, freqmax=15.0, zerophase=False, corners=2, peakf=2.0, wavespeed=2000, sourcelat=None, sourcelon=None, \
                 startt=None, endt=None, units=None ):
        self.stream = Stream()
        self.metric = metric

        try:
            r = [tr.stats.distance for tr in st]
            if not r:
                attach_station_coordinates_from_inventory(inv, st)
                attach_distance_to_stream(st, sourcelat, sourcelon) 
            if not r:
                f"Cannot determine distances from source to stations"
                return
        except:
            f"Cannot determine distances from source to stations"
            return           

        if units=='Counts' and inv: # try to correct          
            for tr in st:
                if tr.stats.channel[2] in 'ENZ' : # filter seismic channels only
                    print('Processing %s' % tr.id)
                    if inv:
                        pre_filt = [freqmin/1,2, freqmin, freqmax, freqmax*1.2]
                        tr.remove_response(output='DISP', inventory=inv, plot=verbose, pre_filt=pre_filt, water_level=60)    
                        self.corrected = True
                        if startt:
                            tr.trim(starttime=startt, endtime=endt)
        elif units=='m':
                print('Units suggest this is already a displacement seismogram. No filtering or instrument correction performed.')
                
        else:
                print('Cannot compute Drs. Need to know units.')
                return

        # We now have Displacement seismogram
        for tr in st:
            this_tr = tr.copy()
            x = (tr.data*100) * np.sqrt(tr.stats.distance*100 * peakf * wavespeed * 100) # everything in cm
            # now we want to reshape the data vector into an array, so we can take advantage of np.mean()
            s = np.size(x) # find the size of the data vector
            nc = int(tr.stats.sampling_rate * sampling_interval) # number of columns
            nr = int(s / nc) # number of rows
            x = x[0:nr*nc] # cut off any trailing samples
            y = x.reshape((nr, nc))
            if verbose:
                print('%s: size %d' % (tr.id, s))
                print('%s: reshaped to %d x %d array (%d samples)' % (tr.id, nr, nc, nr * nc))
            if metric=='mean':
                this_tr.data = np.nanmean(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
            if metric=='median':
                this_tr.data = np.nanmedian(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
            if metric=='max':
                this_tr.data = np.nanmax(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)
            if metric=='rms':
                this_tr.data = np.rms(abs(y),axis=1) # compute mean for each row (this is vectorised; faster than a for loop)     
            this_tr.stats.sampling_rate = 1.0 / sampling_interval # update the sampling rate
            self.stream.append(this_tr)
        self.sampling_interval = sampling_interval



    def plot(self, equal_scale=False, type='linear', percentile=None, linestyle='-'):
        st = self.stream.copy()
        for tr in st:
            tr.data = np.where(tr.data==0, np.nan, tr.data)
        if type=='linear':
            linearplot(st, equal_scale=equal_scale, percentile=percentile, linestyle=linestyle)
        elif type=='log':
            import matplotlib.pyplot as plt
            import matplotlib.dates as mdates
            from math import log2
            plt.rcParams["figure.figsize"] = (10,6)
            fig, ax = plt.subplots()
            for tr in st:
                x = tr.data
                 # now we want to reshape the data vector into an array, so we can take advantage of np.mean()
                s = np.size(x) # find the size of the data vector
                nc = np.max((1, int(log2(s/260))))  # number of columns
                #nc = nc * 4
                nr = int(s / nc) # number of rows
                x = x[0:nr*nc] # cut off any trailing samples
                y = x.reshape((nr, nc))
                y2 = np.nanmax(y,axis=1)    
                t = tr.times("utcdatetime")[::nc]   
                t = [this_t.datetime for this_t in t] 
                #print(t) 
                t = t[:len(y2)]      
                ax.semilogy(t, y2,linestyle, label='%s' % tr.id) #, alpha=0.03)
            ax.format_xdata = mdates.DateFormatter('%H')
            ax.legend()
            plt.xticks(rotation=45)
            plt.ylim((0.2, 100)) # IceWeb plots went from 0.05-30
            plt.yticks([0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0], \
                ['0.2', '0.5', '1', '2', '5', '10', '20', '50', '100'])
            plt.ylabel(r'$D_{RS}$ ($cm^{2}$)')
            plt.xlabel(r'UTC / each point is max $D_{RS}$ in %d minute window' % (tr.stats.delta * nc /60))
            plt.title('Reduced Displacement (%s)\n%s to %s' % (r'$D_{RS}$', t[0].strftime('%d-%b-%Y %H:%M:%S UTC'), t[-1].strftime('%d-%b-%Y %H:%M:%S UTC')))
        return 0   


    def write(self, SDS_TOP):
        DrsSDS_TOP = os.path.join(SDS_TOP,'DRS',self.metric)
        DrsSDSobj = SDS.SDSobj(DrsSDS_TOP, streamobj=self.stream)
        DrsSDSobj.write(overwrite=True)                
    
    def read(self, startt, endt, SDS_TOP, metric='mean', speed=2):
        DrsSDS_TOP = os.path.join(SDS_TOP,'DRS',self.metric)
        thisSDSobj = SDS.SDSobj(DrsSDS_TOP)
        thisSDSobj.read(startt, endt, speed=speed)
        print(thisSDSobj.stream)
        self.stream = thisSDSobj.stream
        self.metric=metric
        self.sampling_interval=self.stream[0].stats.delta  




def linearplot(st, equal_scale=False, percentile=None, linestyle='-'):
    hf = st.plot(handle=True, equal_scale=equal_scale, linestyle=linestyle) #, method='full'); # standard ObsPy plot
    # change the y-axis so it starts at 0
    allAxes = hf.get_axes()
    ylimupper = [ax.get_ylim()[1] for ax in allAxes]
    print(ylimupper)
    if percentile:
        ylimupper = np.array([np.percentile(tr.data, percentile) for tr in st])*1.1
    # if equal_scale True, we set to maximum scale
    print(ylimupper)
    ymax=max(ylimupper)
    for i, ax in enumerate(allAxes):
        if equal_scale==True:
            ax.set_ylim([0, ymax])
        else:
            ax.set_ylim([0, ylimupper[i]])  

def attach_station_coordinates_from_inventory(inventory, st):
    """ attach_station_coordinates_from_inventory """
    for tr in st:
        for netw in inventory.networks:
            for sta in netw.stations:
                if tr.stats.station == sta.code and netw.code == tr.stats.network:
                    for cha in sta.channels:
                        if tr.stats.location == cha.location_code:
                            tr.stats.coordinates = obspy.core.util.AttribDict({
                                'latitude':cha.latitude,
                                'longitude':cha.longitude,
                                'elevation':cha.elevation})
                            
                                                      
def attach_distance_to_stream(st, olat, olon):
    for tr in st:
        try:
            alat = tr.stats['coordinates']['latitude']
            alon = tr.stats['coordinates']['longitude']
            print(alat, alon, olat, olon)
            distdeg = obspy.geodetics.locations2degrees(olat, olon, alat, alon)
            distkm = obspy.geodetics.degrees2kilometers(distdeg)
            tr.stats['distance'] =  distkm * 1000
        except Exception as e:
            print(e)
            print('cannot compute distance for %s' % tr.id)

if __name__ == "__main__":
    pass
