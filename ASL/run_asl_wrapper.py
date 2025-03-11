import os
import glob
import numpy as np
import sys
os.environ["MPLBACKEND"] = "TkAgg"
import matplotlib
matplotlib.use("TkAgg")
from obspy import UTCDateTime, Stream
#from pprint import pprint
script_dir = os.path.dirname(os.path.abspath(__file__))
localLibPath = os.path.join(os.path.dirname(script_dir), 'lib')
sys.path.append(localLibPath)
from libseisGT import detect_network_event
from seisan_classes import set_globals, read_seisandb_apply_custom_function_to_each_event
from SAM import DSAM, DRS
from ASL import ASL, initial_source, make_grid, dome_location
from libDetectionTuner import run_event_detection, signal2noise, plot_detected_stream

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

SEISAN_DATA, DB, station_locationsDF, inv = set_globals(SEISAN_DATA=SEISAN_DATA)

startdate=UTCDateTime(2001,3,1,14,16,0)
enddate=UTCDateTime(2001,3,2)

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
    peakf = kwargs.get("peakf", 8.0)
    metric = kwargs.get("metric", 'VT')
    window_seconds = kwargs.get("window_seconds", 5)
    min_stations = kwargs.get("min_stations", 5)
    outdir =  os.path.join(kwargs.get("outdir", '.'), st[0].stats.starttime.strftime("%Y%m%dT%H%M%S")) 
    interactive = kwargs.get("interactive", False)
    freq = [1.0, 15.0]
    compute_DRS_at_fixed_source = kwargs.get("compute_DRS_at_fixed_source", True)

    if not isinstance(Q, list):
        Q = [Q]
    if not isinstance(peakf, list):        
        peakf = [peakf]
    if not isinstance(surfaceWaveSpeed_kms, list):
        surfaceWaveSpeed_kms = [surfaceWaveSpeed_kms]
    if not isinstance(metric, list):
        metric = [metric]
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
        st.write(os.path.join(outdir, 'stream.mseed'), format='MSEED')

    best_params = {'algorithm': 'zdetect', 'sta': 1.989232720619933, 'lta': 29.97100574322863, 'thr_on': 1.6780449441904004, 'thr_off': 0.6254430984099986}
    if interactive:
        st, best_params, all_params_df = run_event_detection(st, n_trials=500) 
        print(f'best_params: {best_params}')
        all_params_df = all_params_df.sort_values(by='misfit', ascending=True).head(20)
        #all_params_df = all_params_df.nsmallest(20, 'misfit') # just save smallest values
        all_params_df.to_csv(os.path.join(outdir, 'detection_trials.csv'))
        st.write(os.path.join(outdir, 'selected_stream.mseed'), format='MSEED')
        # noise spectra
        snr, fmetrics_dict = signal2noise(st, best_params)
        print(f'tuner signal_to_noise {snr}')
        if fmetrics_dict:
            freq = [fmetrics_dict['f_low'], fmetrics_dict['f_high']]
            print(f'frequency metrics: {fmetrics_dict}')
            peakf = [fmetrics_dict['f_peak']]
        else:
            print('No frequency metrics returned from tuner')

    # detect
    best_trig = detect_network_event(st, minchans=None, threshon=best_params['thr_on'], threshoff=best_params['thr_off'], \
                                     sta=best_params['sta'], lta=best_params['lta'], pad=0.0, best_only=True, freq=freq, \
                                        algorithm=best_params['algorithm'], criterion='cft')
    
    if not best_trig:
        return
    if len(best_trig['trace_ids'])<min_stations:
        return
    print(f'best_trig: {best_trig}')
    
    # subset by detected stations
    detected_st = Stream(traces=[tr for tr in st if tr.id in best_trig['trace_ids']])
    plot_detected_stream(detected_st, best_trig, outfile = os.path.join(outdir, 'detected_stream.png'))

    # noise spectra
    snr, fmetrics_dict = signal2noise(detected_st, best_trig, outfile=os.path.join(outdir, 'signal2noise.png')) # add plot and outfile ?
    print(f'detector signal_to_noise {snr}')

    # compute DSAM data with 1-s time window
    detected_st.trim(starttime=best_trig['time']-1, endtime=best_trig['time']+best_trig['duration']+1)
    dsamObj = DSAM(stream=detected_st, sampling_interval=1.0)
    #print(f'DSAM object for asl_event: {dsamObj}')

    # plot DSAM data
    if dsamObj and isinstance(dsamObj, DSAM):
        dsampng = os.path.join(outdir, 'DSAM.png')
        dsamObj.plot(metrics=metric, equal_scale=True, outfile=dsampng)
        dsamObj.write(outdir, ext='csv')

    # loop over Q, peakf, surfaceWaveSpeed_kms, metric
    print(f'Computing ASL for Q={Q}, peakf={peakf}, surfaceWaveSpeed_kms={surfaceWaveSpeed_kms}, metric={metric}')

    # make grid object
    source = initial_source(lat=dome_location['lat'], lon=dome_location['lon'])
    gridobj = make_grid(center_lat=source['lat'], center_lon=source['lon'], node_spacing_m = 100, grid_size_lat_m = 10000, grid_size_lon_m = 8000)

    for this_Q in Q:
        for this_peakf in peakf:
            for this_surfaceWaveSpeed_kms in surfaceWaveSpeed_kms:
                for this_metric in metric:
                    print(f'Computing ASL for Q={this_Q}, peakf={this_peakf}, surfaceWaveSpeed_kms={this_surfaceWaveSpeed_kms}, metric={this_metric}')

                    if compute_DRS_at_fixed_source:
                        # compute reduced displacement from DSAM data - assumes source at volcano
                        DRSobj = dsamObj.compute_reduced_displacement(inv, source, surfaceWaves=True, Q=this_Q, wavespeed_kms=this_surfaceWaveSpeed_kms, peakf=this_peakf)                
                        # plot DRS data - we assume the frequency band of interest for PDCs is same as the VT band we use for band ratio calculations. I think this is 4-18 Hz. But we could give different banks to the DSAM calculation above
                        DRSmaxrms = DRSobj.max(metric=this_metric)
                        print(f'Maximum DRS assuming fixed source is: {DRSmaxrms}')

                        # plot DRS data
                        if DRSobj and isinstance(DRSobj, DRS):
                            drspng = os.path.join(outdir, f'dome_DRS_Q{this_Q}_f{this_peakf}_v{this_surfaceWaveSpeed_kms}_{this_metric}.png')
                            DRSobj.plot(metrics=this_metric, equal_scale=True, outfile=drspng)               

                    # Create an ASL object with DSAM (displacement amplitude) data, inventory data, and a grid object. The inventory is used for station locations to compute distances
                    aslobj = ASL(dsamObj, this_metric, inv, gridobj, window_seconds)

                    # Compute grid distances
                    aslobj.compute_grid_distances()

                    # Compute amplitude corrections
                    aslobj.compute_amplitude_corrections(surfaceWaves=True, wavespeed_kms=this_surfaceWaveSpeed_kms, Q=this_Q, fix_peakf = this_peakf)

                    # Estimate source location for each time window of window_seconds
                    aslobj.fast_locate()

                    # Print the Event
                    aslobj.print_event()

                    aslobj.save_event(outfile=os.path.join(outdir, f'ASL_Q{this_Q}_f{this_peakf}_v{this_surfaceWaveSpeed_kms}_{this_metric}.qml'))


                    # Plot the source location estimates
                    # originally zoom_level=1, scale=0.2, number=10, equal_size=True
                    aslobj.plot(zoom_level=0, threshold_DR=0.03, \
                                scale=0.2, join=True, number=0, \
                                equal_size=False, add_labels=True, stations=[tr.stats.station for tr in detected_st], \
                                    outfile=os.path.join(outdir, f'ASL_Q{this_Q}_f{this_peakf}_v{this_surfaceWaveSpeed_kms}_{this_metric}.png')) #, show=True)

outdir = os.path.join(os.path.dirname(SEISAN_DATA), 'ASL_DB')
read_seisandb_apply_custom_function_to_each_event(
    startdate, enddate, SEISAN_DATA=SEISAN_DATA, DB='MVOE_', 
    inv=inv, post_process_function=asl_event, verbose=True, 
    bool_clean=True, plot=True, valid_subclasses='re', 
    quality_threshold=1.0, outputType='DISP', freq=[0.5, 30.0], vertical_only=True, 
    max_dropout=4.0,
    # arguments for asl_event follow
    outdir=outdir, Q=23, surfaceWaveSpeed_kms = 1.5, peakf = 8.0, 
    metric='rms', window_seconds=5, min_stations=5,interactive=True,
    )