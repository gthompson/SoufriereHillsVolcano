######################################################################
##   Extended version of ObsPy sdsclient                            ##
######################################################################

from obspy import Stream, read
import obspy.clients.filesystem.sds
import os
import numpy as np
import pandas as pd
import glob


def remove_empty_traces(stream):
    """Removes empty traces, traces full of zeros, and traces full of NaNs from an ObsPy Stream."""
    cleaned_stream = Stream()  # Create a new empty Stream

    for trace in stream:
        # Check if trace is empty (npts == 0)
        if trace.stats.npts == 0:
            continue
        
        # Check if all values are zero
        if np.all(trace.data == 0):
            continue

        # Check if all values are NaN
        if np.all(np.isnan(trace.data)):
            continue

        # If trace passes all checks, add it to the cleaned stream
        cleaned_stream += trace

    return cleaned_stream

# Example usage:
# st = some_obspy_stream
# st = remove_empty_traces(st)


class SDSobj():
    # Call like:
    # myclient = SDS.SDSobj(SDS_TOP, sds_type='D', format='MSEED')
    def __init__(self, SDS_TOP, sds_type='D', format='MSEED', streamobj=Stream()):
        if not os.path.isdir(SDS_TOP):
            os.makedirs(SDS_TOP, exist_ok=True)    
        self.client = obspy.clients.filesystem.sds.Client(SDS_TOP, sds_type=sds_type, format=format)
        self.stream = streamobj
        self.topdir = SDS_TOP

    # Read SDS archive
    def read(self, startt, endt, skip_low_rate_channels=True, trace_ids=None, speed=1, verbose=True ):
        if not trace_ids:
            #nslc_tuples = sdsclient.get_all_nslc(sds_type='D', datetime=startt) 
            trace_ids = self._sds_get_nonempty_traceids(startt, endt)

        st = Stream()
        for trace_id in trace_ids:
            this_st = Stream()
            net, sta, loc, chan = trace_id.split('.')
            if chan[0]=='L' and skip_low_rate_channels:
                print(trace_id,' skipped')
                continue
            if speed==1:
                sdsfiles = self.client._get_filenames(net, sta, loc, chan, startt, endt)
                print(sdsfiles)
                
                for sdsfile in sdsfiles:
                    if os.path.isfile(sdsfile):
                        if verbose:
                            print('st = read("%s")' % sdsfile)
                        that_st = read(sdsfile)
                        print(f'{sdsfile} -> {that_st}')
                        #that_st.merge(method=1,fill_value=0)
                        that_st.merge(method=0) # mark overlaps that disagree with missing samples, and do not fill 
                        #print(sdsfile, that_st)
                        for tr in that_st:
                            this_st.append(tr)  
            elif speed==2:
                this_st = self.client.get_waveforms(net, sta, loc, chan, startt, endt, merge=-1) # use method=0? 
            if this_st:
                print(this_st)
                this_st.trim(startt, endt)   
                #remove_empty_traces(this_st)         
            for tr in this_st:
                st.append(tr)
        if st:
            st.merge(method=0)
            self.stream = st
            return 0
        else:
            self.stream = st
            return 1

        self.stream = st

    # Check SDS archive exists
    def find_which_days_missing(self, stime, etime, net):
        secondsPerDay = 60 * 60 * 24
        dayt = stime
        missing_days = []
        while dayt < etime:
            # check if files for this Julian day already exist
            jday = dayt.strftime('%j')
            yyyy = dayt.strftime('%Y')
            pattern = os.path.join(self.topdir, yyyy, net, '*', '*.D', f"{net}*.{yyyy}.{jday}")
            print('Looking for files in SDS archive matching', pattern)
            existingfiles = glob.glob(pattern)
            print(existingfiles)
            print(f"Found {len(existingfiles)} files")
            if len(existingfiles) == 0:
                missing_days.append(dayt)
            dayt += secondsPerDay
        return missing_days
    
    def _sds_get_nonempty_traceids(self, startday, endday=None, skip_low_rate_channels=True):
        thisday = startday
        trace_ids = []   
        if not endday:
            endday = startday + 86400
        while thisday<endday:
            #print(thisday)
            nslc_tuples = self.client.get_all_nslc(sds_type='D', datetime=thisday) 
            for nslc in nslc_tuples:
                trace_id = '.'.join(nslc)
                if not trace_id in trace_ids:
                    net, sta, loc, chan = nslc
                    if chan[0]=='L' and skip_low_rate_channels:
                        continue
                    if self.client.has_data(net, sta, loc, chan):
                        trace_ids.append(trace_id)  
            thisday = thisday + 86400                    
        return sorted(trace_ids)
 

    def _sds_percent_available_per_day(self, startday, endday, skip_low_rate_channels=True, trace_ids=None, speed=3 ):
        # speed can be 1, 2, or 3. The higher the speed, the less accurate the results.
        #print(trace_ids)
        if not trace_ids:
            trace_ids = self._sds_get_nonempty_traceids(startday, endday)
        if len(trace_ids)==0:
            print('we have no trace_ids')
            return
        else:
            print(trace_ids)
    
        lod = []
        thisday = startday
        print(startday, endday, thisday)
        while thisday<endday:
            print(thisday.date) #, ' ',newline=False)
            thisDict = {}
            thisDict['date']=thisday.date
            for trace_id in trace_ids:
                net, sta, loc, chan = trace_id.split('.')
            
                if speed<3: # quite slow but most accurate
                    sdsfile = self.client._get_filename(net, sta, loc, chan, thisday)
                    this_percent = 0
                    if os.path.isfile(sdsfile):
                        st = obspy.read(sdsfile)
                        st.merge(method=0)
                        tr = st[0]
                        npts_expected = tr.stats.sampling_rate * 86400
                        if speed==1: # slowest, completely accurate
                            npts_got = np.count_nonzero(~np.isnan(tr.data))
                        elif speed==2: # almost as slow, includes np.isnan as valid data sample   
                            npts_got = tr.stats.npts
                        this_percent = 100 * npts_got/npts_expected

                elif speed==3: # much faster
                    this_percent = 100*self.client.get_availability_percentage(net, sta, loc, chan, thisday, thisday+86400)[0]
                
                thisDict[trace_id] = this_percent
            #print(thisDict)
            lod.append(thisDict)
            thisday = thisday + 86400
        availabilityDF = pd.DataFrame(lod)
        print(availabilityDF);
        return availabilityDF, trace_ids

    def plot_availability(self, availabilityDF, outfile=None, FS=12, labels=None):
        import matplotlib.pyplot as plt

        def _get_yticks(s):
            yticklabels = []
            yticks=[]
            for i, el in enumerate(s.to_numpy()):
                yticklabels.append(el)
                yticks.append(i)
            ystep = 1
            if len(yticks)>15:
                ystep=2  
            if len(yticks)>25:
                ystep=3
            if len(yticks)>40:
                ystep=int(len(yticks)*7/100)        
            yticks = yticks[0::ystep]
            yticklabels = yticklabels[0::ystep]
            return (yticks, yticklabels)

        Adf = availabilityDF.iloc[:,1:]/100
        Adata = np.array(Adf, dtype=float) # convert dataframe to numpy array
        if isinstance(labels, list):
            xticklabels = labels
        else:
            xticklabels = Adf.columns
        #print(xticklabels)
        (yticks, yticklabels) = _get_yticks(availabilityDF['date'])
        #print(yticks)

        plt.rc('axes', labelsize=FS) 
        fig = plt.figure(figsize=(FS, FS))
        ax = fig.add_subplot(111)
        #ax.imshow(data, aspect='auto', cmap=plt.cm.gray, interpolation='nearest')
        ax.imshow(np.transpose(1.0-Adata), aspect='auto', cmap=plt.cm.gray, interpolation='nearest')
        #plt.xticks(np.arange(len(xticklabels)), xticklabels)
        plt.yticks(np.arange(len(xticklabels)), xticklabels, fontsize=FS)
        #ax.set_xticklabels(xticklabels, rotation = 90)
        ax.set_yticklabels(xticklabels, rotation = 0, fontsize=FS )
        #plt.yticks(yticks, yticklabels)
        plt.xticks(yticks, yticklabels, rotation = 90, fontsize=FS)
        bottom, top = ax.set_ylim()
        ax.set_ylim(bottom+1, top-1)

        ax.tick_params(axis='both', which='major', labelsize=FS)
        #fig.tight_layout()
        plt.xlabel('Date')
        plt.ylabel('NSLC')
        plt.rc('font', size=FS) 
        plt.rc('xtick', labelsize=FS) 
        plt.rc('ytick', labelsize=FS) 
        plt.grid('on')
        #plt.show()
        if outfile:
            fig.savefig(outfile, dpi=300)
    
    def write(self, overwrite=False, debug=False):
        successful = True
        for tr in self.stream: #RMOVE 0
            if debug:
                if not tr.stats.station == 'MBWH':
                    continue  
                print('/nWant to add this trace:')
                tr.plot();
            sdsfile = self.client._get_filename(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, tr.stats.starttime, 'D')
            sdsdir = os.path.dirname(sdsfile)
            print(sdsdir)       
            if not os.path.isdir(sdsdir):
                os.makedirs(sdsdir, exist_ok=True)
            if overwrite:
                if os.path.isfile(sdsfile): 
                    print('sds.write: ',sdsfile,' already contains data. overwriting')
                else:
                    print('sds.write: ',sdsfile,' does not already exist. Writing')
                try:
                    tr.write(sdsfile, 'mseed')
                except:
                    successful = False
            else:
                if os.path.isfile(sdsfile): # try to load and merge data from file if it already exists
                    st_before = obspy.read(sdsfile)
                    if debug:
                        print('SDS before:')
                        st_before.plot();
                    st_new = st_before.copy()
                    if not st_before == tr:
                        print(f"sds.write: {sdsfile} already contains data: {st_before}")
                        print(f"sds.write: trying to merge new trace {tr}")
                        tr_now = tr.copy()
                        st_new.append(tr_now)
                        try:
                            #st_new.merge(method=0) # overlaps will be marked as gaps, and gaps will not be filled (no fill_value)
                            st_new.merge(method=1, interpolation_samples=0, fill_value=0) # fill_value=0
                            print(f"sds.write: After merge {st_new}")
                        except:
                            print('sds.write: Could not merge. No new data written.')
                            successful = False
                            break
                        if debug:
                            print('SDS new:')
                            st_new.plot();
                        if len(st_new)==1: # merge succeeded
                            try:
                                st_new[0].write(sdsfile, 'mseed')
                            except:
                                successful = False 
                                print('sds.write: after merge, write attempt failed')
                                print(st_new[0])
                                # merge failed, so keep whichever has most samples?
                                if st_before[0].stats.sampling_rate >= 1: 
                                    if tr.stats.sampling_rate >= 1:
                                        if st_before[0].stats.npts < tr.stats.npts:
                                            try:
                                                tr.write(sdsfile, 'mseed')
                                            except:
                                                successful = False
                                                print('sds.write: after merge write fail, attempt to write new trace instead failed')
                                                print(tr)
                        else:
                            successful = False
                            print('sds.write: Cannot write Stream with more than 1 trace to a single SDS file')
                else:    
                    print('sds.write: ',sdsfile,' does not already exist. Writing')
                    try:
                        tr.write(sdsfile, 'mseed')
                    except:
                        successful = False
        return successful
    
    def __str__(self):
        print(f"client={self.client}, stream={self.stream}")


