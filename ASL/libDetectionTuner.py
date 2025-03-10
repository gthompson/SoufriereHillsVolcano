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
    #best_misfit = float("inf")
    min_triggers = 4
    trials_complete = 0
    lod = []
    minchans = len(stream)
    best_misfit = np.Inf
    while trials_complete < n_trials and best_misfit > max_allowed_misfit:
        algorithm = random.choice(algorithms)
        sta = random.uniform(0.1, 5)
        lta = random.uniform(sta * 4, sta * 25)
        thr_on = random.uniform(1.5, 5)
        thr_off = random.uniform(thr_on / 10, thr_on / 2)
        #ratio = random.uniform(1, 2)
        best_trig=None
        if verbose:
            print(f"Trying {algorithm} with STA={sta}, LTA={lta}, THR_ON={thr_on}, THR_OFF={thr_off}")

        # Run coincidence_trigger
        # coincidence_trigger(trigger_type, thr_on, thr_off, stream, thr_coincidence_sum, \
        # trace_ids=None, max_trigger_length=1000000.0, delete_long_trigger=False, trigger_off_extension=0, details=False, event_templates={}, similarity_threshold=0.7, **options)
        #try:
        if True:
            best_trig = detect_network_event(stream, minchans=minchans, threshon=thr_on, threshoff=thr_off, \
                         sta=sta, lta=lta, pad=0.0, best_only=True, verbose=False, freq=None, algorithm=algorithm, criterion='cft')
            '''
            if algorithm == "zdetect":
                lta = None
                triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta)
            elif algorithm == "carlstatrig":
                triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta, lta=lta, ratio=1, quiet=True)
            else:
                triggers = coincidence_trigger(algorithm, thr_on, thr_off, stream, min_triggers, sta=sta, lta=lta)'
            '''
        #except Exception as e:
        #    print(e)
        #    continue

        if not best_trig:   # No triggers found
            if verbose:
                print("No triggers found.")
            continue
        
        '''
        detected_start_times = [UTCDateTime(trig["time"]) for trig in triggers]
        detected_end_times = [UTCDateTime(trig["time"]) + trig["duration"] for trig in triggers]


        if detected_start_times and detected_end_times:
            detected_start = min(detected_start_times)
            detected_end = max(detected_end_times)
            misfit = abs(detected_start - event_start) + abs(detected_end - event_end)

            print(f"Detected event: {detected_start} to {detected_end}, Misfit: {misfit}")

            if misfit < best_misfit:
                best_misfit = misfit
                best_params = {
                    "algorithm": algorithm,
                    "sta": sta,
                    "lta": lta,
                    "thr_on": thr_on,
                    "thr_off": thr_off,
                    "misfit": misfit,
                }'
        '''

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

        ''' other fields of best_trig
            {'cft_peak_wmean': 19.561900329259956,
            'cft_peaks': [19.535644192544272,
                        19.872432918501264,
                        19.622171410201297,
                        19.217352795792998],
            'cft_std_wmean': 5.4565629691954713,
            'cft_stds': [5.292458320417178,
                        5.6565387957966404,
                        5.7582248973698507,
                        5.1190298631982163],
            'coincidence_sum': 4.0,
            'duration': 4.5299999713897705,
            'stations': ['UH3', 'UH2', 'UH1', 'UH4'],
            'time': UTCDateTime(2010, 5, 27, 16, 24, 33, 190000),
            'trace_ids': ['BW.UH3..SHZ', 'BW.UH2..SHZ', 'BW.UH1..SHZ', 'BW.UH4..SHZ']}
        '''

        df = pd.DataFrame(lod)
        best_params = df.loc[df['misfit'].idxmin()].to_dict()
        best_misfit = best_params['misfit']
    #print(df)
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

    df, best_params = run_monte_carlo(app.selected_traces, event_start, event_end, n_trials)
    out_stream = Stream(traces=app.selected_traces)
    plt.close('all')  # Ensures all figures are closed
    root.withdraw()  # Hide the main Tkinter window
    root.quit()  # Quit the Tkinter mainloop
    root.destroy()  # Forcefully close all Tk windows

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

def compute_amplitude_spectra(stream, plot=False, outfile=None):
    if plot:
        plt.figure(figsize=(10, 6))

    for tr in stream:
        # Get time sampling interval (dt) and number of samples (N)
        dt = tr.stats.delta  # Time step
        N = len(tr.data)  # Number of samples

        # Compute FFT
        fft_vals = np.fft.fft(tr.data)
        freqs = np.fft.fftfreq(N, d=dt)  # Frequency axis

        # Compute amplitude spectrum (absolute value of FFT)
        amplitude_spectrum = np.abs(fft_vals)

        # Plot only positive frequencies (since FFT is symmetric)
        positive_freqs = freqs[:N//2]
        positive_amplitudes = amplitude_spectrum[:N//2]
        if plot:
            plt.plot(positive_freqs, positive_amplitudes, label=tr.id)
        

    if plot:
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude")
        plt.title("Amplitude Spectrum of Seismic Signals")
        plt.legend()
        plt.grid()
        plt.xlim(0, 50)  # Adjust frequency range as needed
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
    else:
        return positive_freqs, positive_amplitudes # need to add these to trace objects instead of returning them

def median_spectrum(spectra):
    """
    Computes the median amplitude spectrum.

    Parameters:
        spectra (2D numpy array): Each row is an amplitude spectrum.

    Returns:
        numpy array: Median spectrum.
    """
    return np.median(spectra, axis=0)

def geometric_mean_spectra(spectra):
    """
    Computes the geometric mean of multiple amplitude spectra.

    Parameters:
        spectra (2D numpy array): Each row is an amplitude spectrum.

    Returns:
        numpy array: Geometric mean spectrum.
    """
    spectra = np.array(spectra)

    # Replace zero values to avoid log issues (smallest positive float)
    spectra[spectra == 0] = np.finfo(float).eps

    # Compute geometric mean across spectra (axis=0 means across rows)
    gm_spectrum = np.exp(np.mean(np.log(spectra), axis=0))
    
    return gm_spectrum


def compute_amplitude_ratios(signal_stream, noise_stream, log_scale=False, smooth_window=None, plot=False, outfile=None, verbose=False, average='geometric'):
    if outfile:
        plot = True
    if plot:
        plt.figure(figsize=(10, 6))

    # Create a dictionary for noise traces for quick lookup
    noise_dict = {tr.id: tr for tr in noise_stream}

    spectral_ratios_list = []
    freqs_list = []
    avg_freqs, avg_spectral_ratio = None, None

    for sig_tr in signal_stream:
        trace_id = sig_tr.id
        if trace_id not in noise_dict:
            print(f"Skipping {trace_id}: No matching noise trace.")
            continue

        noise_tr = noise_dict[trace_id]

        # Ensure both traces have the same length
        min_len = min(len(sig_tr.data), len(noise_tr.data))
        sig_data = sig_tr.data[:min_len]
        noise_data = noise_tr.data[:min_len]

        # Get sampling interval (dt) and number of samples (N)
        dt = sig_tr.stats.delta
        N = min_len  # Use the shortest available length

        # Compute FFT for signal and noise
        fft_signal = np.fft.fft(sig_data)
        fft_noise = np.fft.fft(noise_data)
        freqs = np.fft.fftfreq(N, d=dt)  # Frequency axis

        # Compute amplitude spectrum
        amp_signal = np.abs(fft_signal)
        amp_noise = np.abs(fft_noise)

        # Avoid division by zero by replacing zeros with a small number
        amp_noise[amp_noise == 0] = 1e-10  

        # Compute amplitude ratio
        amplitude_ratio = amp_signal / amp_noise

        # Smooth the amplitude ratio using a moving average if requested
        if smooth_window:
            kernel = np.ones(smooth_window) / smooth_window
            amplitude_ratio = np.convolve(amplitude_ratio, kernel, mode="same")

        # Store spectral ratios for computing the overall ratio
        spectral_ratios_list.append(amplitude_ratio[:N//2])
        freqs_list.append(freqs[:N//2])

        # Plot only positive frequencies for individual traces
        if plot:
            if log_scale:
                plt.plot(freqs[:N//2], np.log10(amplitude_ratio[:N//2] + 1), label=trace_id, alpha=0.5, linewidth=1)
            else:
                plt.plot(freqs[:N//2], amplitude_ratio[:N//2], label=trace_id, alpha=0.5, linewidth=1)
    
    if verbose:
        print(f"Processed {len(spectral_ratios_list)} traces    ")

    # Compute the overall spectral ratio by summing all individual ratios
    if spectral_ratios_list:
        if average == 'median':
            avg_spectral_ratio = median_spectrum(spectral_ratios_list)
        elif average == 'geometric':
            avg_spectral_ratio = geometric_mean_spectra(spectral_ratios_list)
        else:
            avg_spectral_ratio = np.mean(np.array(spectral_ratios_list), axis=0)  # Compute the mean spectral ratio
        avg_freqs = np.array(freqs_list[0])  # All traces should have the same frequency bins



        # Plot the overall spectral ratio with a thicker line
        if plot:
            if log_scale:
                plt.plot(avg_freqs, np.log10(avg_spectral_ratio + 1), color="black", linewidth=3, label="Overall Ratio")
            else:
                plt.plot(avg_freqs, avg_spectral_ratio, color="black", linewidth=3, label="Overall Ratio")
    if plot:
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Amplitude Ratio (Signal/Noise)")
        plt.title("Amplitude Ratio of Signal to Noise")
        plt.legend()
        plt.grid()
        plt.xlim(0, 50)  # Adjust frequency range as needed
        plt.ylim(bottom=0)  # Ensure no negative values
        if outfile:
            plt.savefig(outfile)
        else:
            plt.show()
    return avg_freqs, avg_spectral_ratio # could also add to each indivdual trace object as well as returning average

def signal2noise(detected_st, best_trig, outfile=None):
    signal_st = detected_st.copy().trim(starttime=best_trig['time'], endtime=best_trig['time']+best_trig['duration'])
    noise_st = detected_st.copy().trim(endtime=best_trig['time'])
    snr = []
    for i, tr in enumerate(signal_st):
        snr.append(np.nanstd(tr.data) / np.nanstd(noise_st[i].data))
    avg_freqs, avg_spectral_ratio = compute_amplitude_ratios(signal_st, noise_st, log_scale=False, smooth_window=5, outfile=outfile, verbose=True, average='geometric') # add return values, and a plot and outfile 
    if isinstance(avg_freqs, np.ndarray) and isinstance(avg_spectral_ratio, np.ndarray):
        fmetrics_dict = get_bandwidth(avg_freqs, avg_spectral_ratio, threshold=0.5) 
        return snr, fmetrics_dict
    else:
        return snr, None

import numpy as np
import scipy.signal

def get_bandwidth(frequencies, amplitude_ratio, threshold=0.707):
    """
    Estimates the peak frequency, bandwidth, and cutoff frequencies
    from a smoothed amplitude ratio spectrum.

    Parameters:
        frequencies (numpy array): Frequency values (Hz)
        amplitude_ratio (numpy array): Amplitude ratio spectrum (signal/noise)

    Returns:
        dict: Contains 'f_peak', 'f_low', 'f_high', and 'bandwidth'.
    """

    # Step 1: Smooth the amplitude ratio spectrum using a moving average filter
    smoothed_ratio = scipy.signal.savgol_filter(amplitude_ratio, window_length=9, polyorder=2)

    # Step 2: Find the peak frequency (max amplitude in smoothed spectrum)
    peak_index = np.argmax(smoothed_ratio)
    f_peak = frequencies[peak_index]
    A_peak = smoothed_ratio[peak_index]

    # Step 3: Find the -3 dB cutoff points (i.e., where amplitude is 0.707 * A_peak)
    threshold = A_peak * 0.707

    # Find lower frequency cutoff (f_low)
    lower_indices = np.where(smoothed_ratio[:peak_index] < threshold)[0]
    f_low = frequencies[lower_indices[-1]] if len(lower_indices) > 0 else frequencies[0]  # First valid point

    # Find upper frequency cutoff (f_high)
    upper_indices = np.where(smoothed_ratio[peak_index:] < threshold)[0]
    f_high = frequencies[peak_index + upper_indices[0]] if len(upper_indices) > 0 else frequencies[-1]

    # Compute bandwidth
    bandwidth = f_high - f_low

    return {
        "f_peak": f_peak,
        "f_low": f_low,
        "f_high": f_high,
        "bandwidth": bandwidth
    }    

if __name__ == "__main__":
    stream = read('test.mseed', format='MSEED')  # Modify to load actual data
    selected_stream, picked_times, best_params = run_event_detection(stream, n_trials=50)

    print(f"\n Selected {len(selected_stream)} traces")
    print(f" Event Start: {picked_times[0]}, Event End: {picked_times[1]}")
    print(f" Best Parameters: {best_params}")

