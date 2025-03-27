from obspy import Stream, read
import obspy.clients.filesystem.sds
import os
import numpy as np
import pandas as pd
import glob


def remove_empty_traces(stream):
    """
    Removes empty traces, traces full of zeros, and traces full of NaNs from an ObsPy Stream.
    """
    return Stream(tr for tr in stream
                  if tr.stats.npts > 0 and not np.all(tr.data == 0) and not np.all(np.isnan(tr.data)))


class SDSobj:
    """
    FLOVOpy SDS object: Extended SDS client to support reading, writing,
    and plotting SDS archive data with support for trace merging and availability metrics.
    """
    def __init__(self, SDS_TOP, sds_type='D', format='MSEED', streamobj=None):
        os.makedirs(SDS_TOP, exist_ok=True)
        self.client = obspy.clients.filesystem.sds.Client(SDS_TOP, sds_type=sds_type, format=format)
        self.stream = streamobj or Stream()
        self.topdir = SDS_TOP

    def read(self, startt, endt, skip_low_rate_channels=True, trace_ids=None, speed=1, verbose=True):
        if not trace_ids:
            trace_ids = self._sds_get_nonempty_traceids(startt, endt, skip_low_rate_channels)

        st = Stream()
        for trace_id in trace_ids:
            net, sta, loc, chan = trace_id.split('.')
            if chan.startswith('L') and skip_low_rate_channels:
                continue
            try:
                if speed == 1:
                    sdsfiles = self.client._get_filenames(net, sta, loc, chan, startt, endt)
                    for sdsfile in sdsfiles:
                        if os.path.isfile(sdsfile):
                            that_st = read(sdsfile)
                            that_st.merge(method=0)
                            st += that_st
                elif speed == 2:
                    st += self.client.get_waveforms(net, sta, loc, chan, startt, endt, merge=-1)
            except Exception as e:
                if verbose:
                    print(f"Failed to read {trace_id}: {e}")

        st.trim(startt, endt)
        st = remove_empty_traces(st)
        st.merge(method=0)
        self.stream = st
        return 0 if len(st) else 1

    def find_which_days_missing(self, stime, etime, net):
        missing_days = []
        dayt = stime
        while dayt < etime:
            jday = dayt.strftime('%j')
            yyyy = dayt.strftime('%Y')
            pattern = os.path.join(self.topdir, yyyy, net, '*', '*.D', f"{net}*.{yyyy}.{jday}")
            existingfiles = glob.glob(pattern)
            if not existingfiles:
                missing_days.append(dayt)
            dayt += 86400
        return missing_days

    def _sds_get_nonempty_traceids(self, startday, endday=None, skip_low_rate_channels=True):
        endday = endday or startday + 86400
        trace_ids = set()
        thisday = startday
        while thisday < endday:
            for net, sta, loc, chan in self.client.get_all_nslc(sds_type='D', datetime=thisday):
                if chan.startswith('L') and skip_low_rate_channels:
                    continue
                if self.client.has_data(net, sta, loc, chan):
                    trace_ids.add(f"{net}.{sta}.{loc}.{chan}")
            thisday += 86400
        return sorted(trace_ids)

    def _sds_percent_available_per_day(self, startday, endday, skip_low_rate_channels=True, trace_ids=None, speed=3):
        trace_ids = trace_ids or self._sds_get_nonempty_traceids(startday, endday, skip_low_rate_channels)
        lod = []
        thisday = startday
        while thisday < endday:
            row = {'date': thisday.date}
            for trace_id in trace_ids:
                net, sta, loc, chan = trace_id.split('.')
                percent = 0
                try:
                    if speed < 3:
                        sdsfile = self.client._get_filename(net, sta, loc, chan, thisday)
                        if os.path.isfile(sdsfile):
                            st = read(sdsfile)
                            st.merge(method=0)
                            tr = st[0]
                            expected = tr.stats.sampling_rate * 86400
                            npts = np.count_nonzero(~np.isnan(tr.data)) if speed == 1 else tr.stats.npts
                            percent = 100 * npts / expected
                    else:
                        percent = 100 * self.client.get_availability_percentage(net, sta, loc, chan, thisday, thisday + 86400)[0]
                except:
                    percent = 0
                row[trace_id] = percent
            lod.append(row)
            thisday += 86400
        return pd.DataFrame(lod), trace_ids

    def plot_availability(self, availabilityDF, outfile=None, FS=12, labels=None):
        import matplotlib.pyplot as plt
        Adf = availabilityDF.iloc[:, 1:] / 100
        Adata = Adf.to_numpy()
        xticklabels = labels or Adf.columns
        yticks = list(range(len(availabilityDF)))
        yticklabels = availabilityDF['date'].astype(str)

        step = max(1, len(yticks) // 25)
        yticks = yticks[::step]
        yticklabels = yticklabels[::step]

        plt.figure(figsize=(FS, FS))
        ax = plt.gca()
        ax.imshow(1.0 - Adata.T, aspect='auto', cmap='gray', interpolation='nearest')
        ax.set_xticks(yticks)
        ax.set_xticklabels(yticklabels, rotation=90, fontsize=FS)
        ax.set_yticks(np.arange(len(xticklabels)))
        ax.set_yticklabels(xticklabels, fontsize=FS)
        ax.set_xlabel('Date')
        ax.set_ylabel('NSLC')
        ax.grid(True)
        if outfile:
            plt.savefig(outfile, dpi=300)

    def write(self, overwrite=False, debug=False):
        successful = True
        for tr in self.stream:
            try:
                sdsfile = self.client._get_filename(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, tr.stats.starttime, 'D')
                os.makedirs(os.path.dirname(sdsfile), exist_ok=True)

                if not overwrite and os.path.isfile(sdsfile):
                    existing = read(sdsfile)
                    new_st = existing.copy().append(tr).merge(method=1, fill_value=0)
                    if len(new_st) == 1:
                        new_st[0].write(sdsfile, format='MSEED')
                    else:
                        raise ValueError("Cannot write Stream with more than 1 trace to a single SDS file")
                else:
                    tr.write(sdsfile, format='MSEED')

            except Exception as e:
                print(f"Write failed for {tr.id}: {e}")
                successful = False

        return successful

    def __str__(self):
        return f"client={self.client}, stream={self.stream}"
