import os
import glob
import numpy as np
import matplotlib
import sys

# Use 'MacOSX' backend if on macOS
if sys.platform == "darwin":
    matplotlib.use("MacOSX")  #  Works best for macOS
    #matplotlib.use("Qt5Agg")

# If 'MacOSX' still crashes, use 'Qt5Agg'
elif sys.platform.startswith("linux"):
    matplotlib.use("Qt5Agg")  #  Linux safe choice

# Now import Matplotlib
import matplotlib.pyplot as plt
from obspy import read_inventory, UTCDateTime
from pprint import pprint

localLibPath = os.path.join('..', 'lib')
sys.path.append(localLibPath)
from libseisGT import detect_network_event
from seisan_classes import set_globals, read_seisandb_apply_custom_function_to_each_event
#import libSeisan2Pandas as seisan
#import libMVO 
from SAM import DSAM
from ASL import ASL, initial_source, make_grid
from InventoryTools import show_response #, has_response
from obspy import Stream

#SDS_DIR = '/data/SDS' #'/Volumes/NTFS_2TB_EXT/SDS'
HOME = os.path.expanduser('~')
if sys.platform == "darwin":
    print("Running on macOS")
    if os.path.exists('/Volumes/DATA'):
        SEISAN_DATA = os.path.join('/Volumes/DATA/SEISAN_DB')
    else:
        raise Exception('No SEISAN_DB found: did you forget to mount the DATA drive?')
elif sys.platform.startswith("linux"):
    print("Running on Linux")
    SEISAN_DATA = '/data/SEISAN_DB'
#DATA_DIR = os.path.join(HOME, 'Dropbox/BRIEFCASE/MESS2024/skience2024_GTplus/02 Volcano Monitoring/data')
#SAM_DIR = os.path.join(DATA_DIR,'continuous','SAM')
#SAM_DIR = '/data/SAM' #os.path.join(DATA_DIR,'continuous','SAM')
#DEM_DIR = os.path.join(DATA_DIR,'DEMs')
#RESPONSE_DIR = os.path.join(DATA_DIR,'responses')
#SEISAN_DATA = os.path.join( '/data', 'SEISAN_DB')
#DB = 'MVOE_'

SEISAN_DATA, DB, station_locationsDF, inv = set_globals(SEISAN_DATA=SEISAN_DATA)
print(inv)
show_response(inv)

startdate=UTCDateTime(2001,3,1,14,16,0)
enddate=UTCDateTime(2001,3,2)
import numpy as np
import matplotlib.pyplot as plt

def plot_amplitude_spectra(stream):
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

        plt.plot(positive_freqs, positive_amplitudes, label=tr.id)

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude")
    plt.title("Amplitude Spectrum of Seismic Signals")
    plt.legend()
    plt.grid()
    plt.xlim(0, 50)  # Adjust frequency range as needed
    plt.show()




def plot_amplitude_ratios(signal_stream, noise_stream, log_scale=False, smooth_window=None):
    plt.figure(figsize=(10, 6))

    # Create a dictionary for noise traces for quick lookup
    noise_dict = {tr.id: tr for tr in noise_stream}

    spectral_ratios_list = []
    freqs_list = []

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
        if log_scale:
            plt.plot(freqs[:N//2], np.log10(amplitude_ratio[:N//2] + 1), label=trace_id, alpha=0.5, linewidth=1)
        else:
            plt.plot(freqs[:N//2], amplitude_ratio[:N//2], label=trace_id, alpha=0.5, linewidth=1)

    # Compute the overall spectral ratio by summing all individual ratios
    if spectral_ratios_list:
        avg_spectral_ratio = np.mean(np.array(spectral_ratios_list), axis=0)  # Compute the mean spectral ratio
        avg_freqs = freqs_list[0]  # All traces should have the same frequency bins

        # Plot the overall spectral ratio with a thicker line
        if log_scale:
            plt.plot(avg_freqs, np.log10(avg_spectral_ratio + 1), color="black", linewidth=3, label="Overall Ratio")
        else:
            plt.plot(avg_freqs, avg_spectral_ratio, color="black", linewidth=3, label="Overall Ratio")

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Amplitude Ratio (Signal/Noise)")
    plt.title("Amplitude Ratio of Signal to Noise")
    plt.legend()
    plt.grid()
    plt.xlim(0, 50)  # Adjust frequency range as needed
    plt.ylim(bottom=0)  # Ensure no negative values
    plt.show()

from libDetectionTuner import run_event_detection
def asl_event(st, raw_st, **kwargs):
    print(f"Received kwargs in asl_event: {kwargs}")
    if len(st) > 0 and isinstance(st, Stream):
        pass
    else:
        print(f"Empty stream for asl_event: {st}")
        return

    # Set default values for parameters (if not provided in **kwargs)
    Q = kwargs.get("Q", 23)  
    surfaceWaveSpeed_kms = kwargs.get("surfaceWaveSpeed_kms", 2.5)
    peakf = kwargs.get("peakf", 10.0)
    metric = kwargs.get("metric", 'VT')
    window_seconds = kwargs.get("window_seconds", 5)
    min_stations = kwargs.get("min_stations", 5)
    outdir =  os.path.join(kwargs.get("outdir", '.'), st[0].stats.starttime.strftime("%Y%m%dT%H%M%S")) 
    interactive = kwargs.get("interactive", False)
    os.makedirs(outdir, exist_ok=True)

    print(f"Using Q={Q}, surfaceWaveSpeed_kms={surfaceWaveSpeed_kms}, peakf={peakf}")    
    if raw_st and isinstance(raw_st, Stream):
        rawstreampng = os.path.join(outdir, 'rawstream.png')
        raw_st.plot(equal_scale=False, outfile=rawstreampng);

    # plot pre-processed stream
    print(f'Stream for asl_event: {st}')
    if st and isinstance(st, Stream):
        streampng = os.path.join(outdir, 'stream.png')
        st.plot(equal_scale=False, outfile=streampng);
    
    # plot amplitude spectra
    #plot_amplitude_spectra(st)    

    # --- Run the pipeline ---

    #st.write('test.mseed', format='MSEED')
    best_params = {'algorithm': 'zdetect', 'sta': 1.989232720619933, 'lta': 29.97100574322863, 'thr_on': 1.6780449441904004, 'thr_off': 0.6254430984099986}
    if interactive:
        st, picked_times, best_params = run_event_detection(st, n_trials=50) 
        print(st)

    # detect
    best_trig = detect_network_event(st, minchans=None, threshon=best_params['thr_on'], threshoff=best_params['thr_off'], \
                                     sta=best_params['sta'], lta=best_params['lta'], pad=0.0, best_only=True, freq=[2.0, 15.0], \
                                        algorithm=best_params['algorithm'])
    if not best_trig:
        return
    pprint(best_trig)
    if len(best_trig['trace_ids'])<min_stations:
        return
    
    # subset by detected stations
    detected_st = Stream(traces=[tr for tr in st if tr.id in best_trig['trace_ids']])

    # plot detected stream
    trig_time = best_trig['time'].matplotlib_date
    trig_end = (best_trig['time'] + best_trig['duration']).matplotlib_date
    fig = detected_st.plot(show=False)  # show=False prevents it from auto-displaying
    for ax in fig.axes:  # Stream.plot() returns multiple axes (one per trace)
        ax.axvline(trig_time, color='r', linestyle='--', label="Trigger Start")
        ax.axvline(trig_end, color='b', linestyle='--', label="Trigger End")
        ax.legend()
    #plt.savefig() # Display the plot    
    fig.savefig(os.path.join(outdir, 'detectedstream.png'), dpi=300, bbox_inches="tight")

    # noise spectra
    signal_st = detected_st.copy().trim(starttime=best_trig['time'], endtime=best_trig['time']+best_trig['duration'])
    noise_st = detected_st.copy().trim(endtime=best_trig['time'])
    plot_amplitude_spectra(noise_st) 
    plot_amplitude_ratios(signal_st, noise_st, log_scale=False, smooth_window=5)
    

    # compute DSAM data with 1-s time window
    dsamObj = DSAM(stream=detected_st, sampling_interval=1.0)
    #print(f'DSAM object for asl_event: {dsamObj}')

    # plot DSAM data
    if dsamObj and isinstance(dsamObj, DSAM):
        dsampng = os.path.join(outdir, 'DSAM.png')
        dsamObj.plot(metrics=metric, equal_scale=True, outfile=dsampng)

    '''
    # compute reduced displacement from DSAM data - assumes source at volcano
    DRSobj = dsamObj.compute_reduced_displacement(inv, source, surfaceWaves=True, Q=Q, wavespeed_kms=surfaceWaveSpeed_kms, peakf=peakf)

    # plot DRS data - we assume the frequency band of interest for PDCs is same as the VT band we use for band ratio calculations. I think this is 4-18 Hz. But we could give different banks to the DSAM calculation above
    DRSmaxrms = DRSobj.max(metric=metric)
    print(f'Maximum DRS assuming fixed source is: {DRSmaxrms}')

    # plot DRS data
    if DRSobj and isinstance(DRSobj, DRS):
        drspng = os.path.join(outdir, 'DRS.png')
        DRSobj.plot(metrics=metric, equal_scale=True, outfile=drspng)
    '''

    # same grid as before
    source = initial_source()
    gridobj = make_grid()

    # Create an ASL object with DSAM (displacement amplitude) data, inventory data, and a grid object. The inventory is used for station locations to compute distances
    aslobj = ASL(dsamObj, 'VT', inv, gridobj, window_seconds)

    # Compute grid distances
    aslobj.compute_grid_distances()

    # Compute amplitude corrections
    aslobj.compute_amplitude_corrections(surfaceWaves=True, wavespeed_kms=surfaceWaveSpeed_kms, Q=Q, fix_peakf = peakf)

    # Estimate source location for each time window of window_seconds
    source_pf = aslobj.fast_locate()
    #source_pf = aslobj.locate()
    #print(source_pf)

    # Plot the source location estimates
    # originally zoom_level=1, scale=0.2, number=10, equal_size=True
    aslobj.plot(source_pf, zoom_level=0, threshold_DR=0.03, \
                scale=0.2, join=True, number=0, \
                equal_size=False, add_labels=True, outfile=os.path.join(outdir, 'ASL.png')) 
                

read_seisandb_apply_custom_function_to_each_event(startdate, enddate, 
                                                SEISAN_DATA=SEISAN_DATA, 
                                                DB='MVOE_', 
                                                inv=inv, 
                                                post_process_function=asl_event, 
                                                #post_process_function=None, 
                                                verbose=True, 
                                                bool_clean=True, 
                                                plot=True, 
                                                valid_subclasses='re', 
                                                quality_threshold=1.0, 
                                                outputType='DISP', 
                                                freq=[0.5, 30.0],
                                                vertical_only=True, 
                                                # arguments for asl_event follow
                                                outdir='.',
                                                Q=23, 
                                                surfaceWaveSpeed_kms = 1.5,
                                                peakf = 8.0, 
                                                metric='rms',
                                                window_seconds=5, 
                                                min_stations=5,
                                                interactive=True,
                                                )
                   