import numpy as np
import tkinter as tk
from tkinter import messagebox
from obspy import read, UTCDateTime, Stream
#from obspy.signal.trigger import coincidence_trigger
import random
import os
os.environ["MPLBACKEND"] = "TkAgg"
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from pprint import pprint
import pandas as pd
from libseisGT import detect_network_event
from metrics import estimate_snr, compute_amplitude_spectra, plot_amplitude_spectra, \
    compute_amplitude_ratios, plot_amplitude_ratios, get_bandwidth

class SeismicGUI:
    def __init__(self, master, stream):
        self.master = master
        self.master.title("Seismic Event Detection")

        self.stream = stream
        self.selected_traces = stream.copy()  # Default: All traces selected
        self.picked_times = []
        self.mode = "select_traces"

        self.frame_left = tk.Frame(master)
        self.frame_left.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)

        self.vars = []
        self.checkboxes = []
        for i, tr in enumerate(stream):
            var = tk.BooleanVar(value=True)  # Default to checked
            chk = tk.Checkbutton(self.frame_left, text=tr.id, variable=var)
            chk.pack(anchor="w")
            self.checkboxes.append(chk)
            self.vars.append(var)

        self.fig, self.axs = plt.subplots(len(stream), 1, figsize=(10, 6), sharex=True)
        plt.subplots_adjust(hspace=0.3)

        self.canvas = FigureCanvasTkAgg(self.fig, master)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.plot_traces()

        self.button = tk.Button(self.frame_left, text="Select Traces", command=self.select_traces)
        self.button.pack(pady=10)

        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('button_press_event', self.on_mouse_click)

    def plot_traces(self):
        """Plots traces separately in linked subplots, ensuring correct time scaling and linking."""
        for ax in self.axs:
            ax.clear()

        start_time = min(tr.stats.starttime for tr in self.stream)
        end_time = max(tr.stats.endtime for tr in self.stream)

        start_time_num = mdates.date2num(start_time.datetime)
        end_time_num = mdates.date2num(end_time.datetime)

        if self.mode == "select_traces":
            self.fig.suptitle("Select traces using checkboxes, then click 'Select Traces'")

            for i, tr in enumerate(self.stream):
                absolute_times = mdates.date2num([start_time + t for t in tr.times()])
                data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
                self.axs[i].plot(absolute_times, data_norm, label=tr.id)
                self.axs[i].legend(loc="upper right")
                self.axs[i].set_ylabel(tr.id)

            self.axs[-1].set_xlabel("Time (UTC)")

        elif self.mode == "pick_times":
            self.fig.suptitle("Click to select event start and end times")

            for i, tr in enumerate(self.selected_traces):
                absolute_times = mdates.date2num([tr.stats.starttime + t for t in tr.times()])
                data_norm = tr.data / max(abs(tr.data)) if max(abs(tr.data)) != 0 else tr.data
                self.axs[i].plot(absolute_times, data_norm, label=tr.id)
                self.axs[i].legend(loc="upper right")
                self.axs[i].set_ylabel(tr.id)

            self.axs[-1].set_xlabel("Time (UTC)")

            self.axs[0].set_xlim(start_time_num, end_time_num)

            self.cursor_lines = [ax.axvline(x=start_time_num, color='r', linestyle='dotted', lw=1) for ax in self.axs]

        for ax in self.axs:
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M:%S"))

        self.canvas.draw()

    def select_traces(self):
        """Handles trace selection and switches to event picking mode."""
        self.selected_traces = [self.stream[i] for i, var in enumerate(self.vars) if var.get()]

        if not self.selected_traces:
            messagebox.showerror("Error", "Please select at least one trace!")
            return

        self.mode = "pick_times"
        self.button.config(text="Pick Start and End Times", state=tk.DISABLED)
        self.plot_traces()

    def on_mouse_move(self, event):
        """Move dotted cursor in pick_times mode across all subplots."""
        if self.mode == "pick_times" and event.xdata:
            for cursor_line in self.cursor_lines:
                if not isinstance(event.xdata, list): # is not None:  # Ensure valid xdata
                    cursor_line.set_xdata([event.xdata]) 
                else:
                    cursor_line.set_xdata(event.xdata)
            self.canvas.draw()

    def on_mouse_click(self, event):
        """Handles event start and end time selection."""
        if self.mode == "pick_times" and event.xdata:
            picked_time = UTCDateTime(mdates.num2date(event.xdata))
            for ax in self.axs:
                ax.axvline(x=event.xdata, color='r', linestyle='solid', lw=2)

            self.picked_times.append(picked_time)
            print(f" Picked event time: {picked_time}")

            if len(self.picked_times) == 1:
                self.fig.suptitle("Click to select event end time")
            elif len(self.picked_times) == 2:
                plt.close(self.fig)  #  Closes only the current figure
                plt.close('all')  #  Ensures all figures are closed                
                self.master.quit()

            self.canvas.draw()

def run_monte_carlo(stream, event_start, event_end, n_trials=5000, verbose=False, max_allowed_misfit=4.0):
    """Runs Monte Carlo simulations with different trigger parameters."""
    print('Finding best autodetect parameters by running Monte Carlo simulations...')
    algorithms = ["recstalta", "classicstalta", "zdetect", "carlstatrig", "delayedstalta"]
    # classicstalta: Classic STA/LTA algorithm
        # classic_sta_lta(tr.data, nsta, nlta) nsta is number of samples
    # recstalta: Recursive STA/LTA algorithm
        # recursive_sta_lta(tr.data, nsta, nlta)
    # zdetect: Z-detect algorithm
        # z_detect(tr.data, nsta) https://docs.obspy.org/packages/autogen/obspy.signal.trigger.z_detect.html#obspy.signal.trigger.z_detect
    # carlstatrig: Carl-Sta-Trig algorithm
        # carl_sta_trig(tr.data, nsta, nlta, ratio, quiet) https://docs.obspy.org/packages/autogen/obspy.signal.trigger.carl_sta_trig.html#obspy.signal.trigger.carl_sta_trig
    # delayed_sta_lta: Delayed STA/LTA algorithm
        # delayed_sta_lta(tr.data, nsta, nlta)
    best_params = None
    trials_complete = 0
    lod = []
    minchans = len(stream)
    best_misfit = np.Inf
    stream_duration = stream[0].stats.endtime - stream[0].stats.starttime
    while trials_complete < n_trials and best_misfit > max_allowed_misfit:
        algorithm = random.choice(algorithms)
        sta = random.uniform(0.1, 5)
        if sta > stream_duration/20:
            sta = stream_duration/20       
        lta = random.uniform(sta * 4, sta * 25)
        if lta > stream_duration/2:
            lta = stream_duration/2
        thr_on = random.uniform(1.5, 5)
        thr_off = random.uniform(thr_on / 10, thr_on / 2)
        #ratio = random.uniform(1, 2)
        best_trig=None
        if verbose:
            print(f"Trying {algorithm} with STA={sta}, LTA={lta}, THR_ON={thr_on}, THR_OFF={thr_off}")

        # Run coincidence_trigger
        if True:
            best_trig = detect_network_event(stream, minchans=minchans, threshon=thr_on, threshoff=thr_off, \
                         sta=sta, lta=lta, pad=0.0, best_only=True, verbose=False, freq=None, algorithm=algorithm, criterion='cft')

        if not best_trig:   # No triggers found
            if verbose:
                print("No triggers found.")
            continue
        
        trials_complete += 1
        best_trig["algorithm"] = algorithm
        best_trig["sta"] = sta
        best_trig["lta"] = lta
        best_trig["thr_on"] = thr_on
        best_trig["thr_off"] = thr_off
        best_trig['endtime'] = best_trig["time"] + best_trig["duration"]
        best_trig['misfit'] = abs(best_trig['time'] - event_start) + abs(best_trig['endtime'] - event_end)
        best_trig['event_start'] = event_start
        best_trig['event_end'] = event_end
        lod.append(best_trig)
        if trials_complete % 10 == 0:
            df = pd.DataFrame(lod)
            best_params = df.loc[df['misfit'].idxmin()].to_dict()
            best_misfit = best_params['misfit']
            print(f'{trials_complete} of {n_trials} trials complete, lowest misfit is {best_misfit}')

    df = pd.DataFrame(lod)
    best_params = df.loc[df['misfit'].idxmin()].to_dict()
    best_misfit = best_params['misfit']
  
    return df, best_params

def run_event_detection(stream, n_trials=50):
    """Runs GUI and returns selected traces, event times, and best parameters."""
    root = tk.Tk()
    app = SeismicGUI(root, stream)
    root.mainloop()

    if len(app.picked_times) != 2:
        print("⚠️ Error: No valid event times selected.")
        return None, None, None, None

    event_start, event_end = app.picked_times
    out_stream = Stream(traces=app.selected_traces)
    plt.close('all')  # Ensures all figures are closed
    root.withdraw()  # Hide the main Tkinter window
    root.quit()  # Quit the Tkinter mainloop
    root.destroy()  # Forcefully close all Tk windows    
    df, best_params = run_monte_carlo(out_stream, event_start, event_end, n_trials)


    return out_stream, best_params, df

def plot_detected_stream(detected_st, best_trig, outfile=None):
    """Plots the detected stream with trigger start and end times."""

    # plot detected stream
    trig_time = best_trig['time'].matplotlib_date
    trig_end = (best_trig['time'] + best_trig['duration']).matplotlib_date
    fig = detected_st.plot(show=False)  # show=False prevents it from auto-displaying
    for ax in fig.axes:  # Stream.plot() returns multiple axes (one per trace)
        ax.axvline(trig_time, color='r', linestyle='--', label="Trigger Start")
        ax.axvline(trig_end, color='b', linestyle='--', label="Trigger End")
        ax.legend()
    if outfile:  
        fig.savefig(outfile, dpi=300, bbox_inches="tight")
    else:
        plt.show()


def detection_snr(detected_st, best_trig, outfile=None, method='std',
                 spectral_smooth=5, spectral_average='geometric',
                 spectral_threshold=0.5, freq_band=None):
    """
    Estimate time-domain and spectral-domain SNR from a triggered waveform stream.

    Parameters
    ----------
    detected_st : obspy.Stream
        The full stream of waveform data containing signal and noise.
    best_trig : dict
        Dictionary with keys 'time' (UTCDateTime) and 'duration' (float).
    outfile : str, optional
        If provided, a spectral amplitude ratio plot will be saved to this file.
    method : str, optional
        Time-domain SNR method: 'std', 'max', 'rms', or 'spectral'.
    spectral_smooth : int, optional
        Moving average window for smoothing amplitude ratios.
    spectral_average : str, optional
        Spectral averaging method: 'geometric', 'mean', 'median'.
    spectral_threshold : float, optional
        Bandwidth threshold as fraction of peak for get_bandwidth().
    freq_band : tuple(float, float), optional
        Frequency band (Hz) to average spectral SNR over.

    Returns
    -------
    snr_list : list of float
        SNR values computed using the specified method.
    fmetrics_dict : dict or None
        Bandwidth metrics from the average spectral ratio, or None if unavailable.
    """
    trig_stime = best_trig['time']
    trig_etime = trig_stime + best_trig['duration']

    # Split into signal and noise windows
    signal_st = detected_st.copy().trim(starttime=trig_stime, endtime=trig_etime)
    noise_st = detected_st.copy().trim(endtime=trig_stime)

    # Estimate time-domain SNR
    snr_list, _, _ = estimate_snr(
        detected_st, method=method, split_time=(trig_stime, trig_etime),
        spectral_kwargs={
            'smooth_window': spectral_smooth,
            'average': spectral_average
        },
        freq_band=freq_band
    )

    # Estimate spectral amplitude ratios
    avg_freqs, avg_spectral_ratio, indiv_ratios, freqs_list, trace_ids, _ = compute_amplitude_ratios(
        signal_st, noise_st,
        smooth_window=spectral_smooth,
        average=spectral_average,
        verbose=True
    )

    # Optional plot
    if isinstance(avg_freqs, np.ndarray) and isinstance(avg_spectral_ratio, np.ndarray):
        if outfile:
            plot_amplitude_ratios(
                avg_freqs, avg_spectral_ratio,
                individual_ratios=indiv_ratios,
                freqs_list=freqs_list,
                trace_ids=trace_ids,
                log_scale=True,
                outfile=outfile,
                threshold=spectral_threshold
            )

        # Bandwidth metrics
        fmetrics_dict = get_bandwidth(avg_freqs, avg_spectral_ratio, threshold=spectral_threshold)
        return snr_list, fmetrics_dict

    else:
        return snr_list, None
  
def real_time_optimization(band='all'):
    corners = 2
    if band=='VT':
        # VT + false
        sta_secs = 1.4
        lta_secs = 7.0
        threshON = 2.4
        threshOFF = 1.2
        freqmin = 3.0
        freqmax = 18.0
    elif band=='LP':
        # LP + false
        sta_secs = 2.3
        lta_secs = 11.5
        threshON = 2.4
        threshOFF = 1.2
        freqmin = 0.8
        freqmax = 10.0        
    elif band=='all':
        # all = LP + VT + false
        sta_secs = 2.3
        lta_secs = 11.5
        threshON = 2.4
        threshOFF = 1.2
        freqmin = 1.5
        freqmax = 12.0        
    threshOFF = threshOFF / threshON
        
    return sta_secs, lta_secs, threshON, threshOFF, freqmin, freqmax, corners

if __name__ == "__main__":
    stream = read('test.mseed', format='MSEED')  # Modify to load actual data
    selected_stream, picked_times, best_params = run_event_detection(stream, n_trials=50)

    print(f"\n Selected {len(selected_stream)} traces")
    print(f" Event Start: {picked_times[0]}, Event End: {picked_times[1]}")
    print(f" Best Parameters: {best_params}")

