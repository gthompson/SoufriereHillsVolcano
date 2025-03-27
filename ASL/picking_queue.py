import os
#import glob
#import numpy as np
import sys
#from obspy import UTCDateTime, Stream
script_dir = os.path.dirname(os.path.abspath(__file__))
localLibPath = os.path.join(os.path.dirname(script_dir), 'lib')
sys.path.append(localLibPath)
from seisan_classes import set_globals, \
    read_seisandb_apply_custom_function_to_each_event
from ASL import dome_location
#from libDetectionTuner import run_event_detection, signal2noise, plot_detected_stream, run_monte_carlo
from metrics import ampengfftmag

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

#startdate=UTCDateTime(2001,3,1,14,16,0)
#enddate=UTCDateTime(2001,3,2)
startdate=UTCDateTime(2001,1,1)
enddate=UTCDateTime(2001,1,2)

def summarize_event(st, raw_st, **kwargs):
    print(f"Received kwargs in asl_event: {kwargs}")
    if len(st) > 0 and isinstance(st, Stream):
        pass
    else:
        print(f"Empty stream for asl_event: {st}")
        return

    # Set default values for parameters (if not provided in **kwargs)
    inventory = kwargs.get('inventory', None)
    outdir =  os.path.join(kwargs.get("outdir", '.'), st[0].stats.starttime.strftime("%Y%m%dT%H%M%S")) 

    os.makedirs(outdir, exist_ok=True)

    # import initial_source from ASL
    source = initial_source(lat=dome_location['lat'], lon=dome_location['lon'])

    # import ampengfftmag from metrics
    newst, df = ampengfftmag(st, inventory, source,
                 model='body', Q=50, c_earth=2500, correction=3.7,
                 a=1.6, b=-0.15, g=0,
                 threshold=0.707, window_length=9, polyorder=2,
                 differentiate=True, verbose=True, snr_method='std',
                 snr_split_time=None, snr_window_length=1.0,
                 snr_min=None)

    """ what we want to do here is:
    2. try to autoclassify event"
    4. add this information to a picking queue dataframe that we modify inside this function?
    """


outdir = os.path.join(os.path.dirname(SEISAN_DATA), 'ASL_DB')
read_seisandb_apply_custom_function_to_each_event(
    startdate, enddate, SEISAN_DATA=SEISAN_DATA, DB='MVOE_', 
    inv=inv, post_process_function=summarize_event, verbose=True, 
    bool_clean=True, plot=True, valid_subclasses='re', 
    quality_threshold=1.0, outputType='DISP', freq=[0.2, 25.0], vertical_only=True, 
    max_dropout=4.0,
    # arguments for summarize_event follow
    outdir=outdir,
    )

"""
to do:

1. make plots from the dataframe
    * magnitude vs time

"""