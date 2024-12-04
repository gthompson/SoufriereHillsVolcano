import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

##### imports for functions from 510_Instrumental_VEI
#import os
import sys
import obspy
sys.path.append('lib')
from SAM import DSAM, VSEM, DR, DRS

# functions already here

from math import log10, floor
def sigfigs(x, n=2):
    if x==0:
        return 0
    else:
        return round(x, n-int(floor(log10(abs(x))))-1)

def summarize_catdf(cat):
    cat = cat.loc[:, ~cat.columns.str.contains('^Unnamed')]
    display(cat)
    print(f'# events = {len(cat)}' + '\n\n')
    print(f'columns = {cat.columns}')

def linregress(catdf, col1, col2, plot=False, print_stats=False, minrows=5, r2min=0.3, show_now=True):

    for col in [col1, col2]:
        if not col in catdf.columns:
            print(f'linregress: {col} not in {catdf.columns}')
            return 0,0,0
    
    model = LinearRegression()

    df = catdf.copy()
    df = df[[col1, col2]]
    df[col1]=df[col1].astype(float)
    df[col2]=df[col2].astype(float)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    #df.replace(NaN, np.nan, inplace=True)
    
    df.dropna(inplace=True)  
    
    for i, row in df.iterrows():
        if np.isnan(row[col2]):
            df.drop(i, inplace=True)
    
    #display(df)
    
    if len(df)<minrows:
        print(f'{col2} vs. {col1}: Not enough rows')
        return 0,0,0,0,0
        


    x = df[col1].to_numpy().reshape(-1, 1)
    y = df[col2].to_numpy().reshape(-1, 1)

    # Train the model
    if len(x)==0 or len(y)==0:
        return 0,0,0,0,0
    else:
        print('\n', f'Got {len(df)} events')
        model.fit(x, y)
    #print(len(x), len(y))

    # Evaluate the model
    r2 = model.score(x, y)
    y_pred = model.predict(x)
    c = model.intercept_[0]
    m = model.coef_[0][0]
    #print(f'{col2} = {c} + {m} {col1}')
    # y = mx + c implies x = (y-c)/m = y/m - c/m
    if c >= 0.0:
        print(f"{col2} = {sigfigs(m,n=3)} {col1} + {sigfigs(c,n=3)}")
        print(f"alt: {col1} = {sigfigs(1/m,n=3)} {col2} - {sigfigs(c/m,n=3)}")       
    else:
        print(f"{col2} = {sigfigs(m,n=3)} {col1} - {sigfigs(c*-1,n=3)}")   
        print(f"alt: {col1} = {sigfigs(1/m,n=3)} {col2} + {sigfigs(c*-1/m,n=3)}") 
        
    if print_stats:
        #print(f"R-squared value: {r2}")
        print(f"R-squared value: {sigfigs(r2,n=3)}")
        # The mean squared error
        print("Mean squared error: %.2f" % mean_squared_error(y, y_pred))
        # The coefficient of determination: 1 is perfect prediction
        #print("Coefficient of determination: %.2f" % r2_score(y, y_pred))

    # Plot outputs
    if plot and r2>r2min:
        fig, ax = plt.subplots()
        plt.scatter(x, y, color="black")
        plt.plot(x, y_pred, color="blue", linewidth=3)
        plt.xlabel(col1)
        plt.ylabel(col2)
        if show_now:
            plt.show()
        
    else:
        fig=None
        ax=None
    return [sigfigs(model.intercept_[0],n=3), sigfigs(model.coef_[0][0],n=3), r2, fig, ax]

def linearRegressMagnitudesBySubclass(cat, subclass_col='subclass', mag_columns=['ML', 'ME'], \
                                      subclasses=['r', 'e', 'l', 'h', 't'], \
                                      plot=False, print_stats=False, mfixed=None):



    print(f'# events = {len(cat)}')
    for i,col1 in enumerate(mag_columns):
        for i2,col2 in enumerate(mag_columns): # alt is for col2 in mag_columns[i+1:]:
            if not col1==col2:
                if mfixed:
                    fix_slope(cat, col1, col2, mfixed=mfixed[i2]/mfixed[i], plot=plot, print_stats=print_stats)
                else:
                    linregress(cat, col1, col2, plot=plot, print_stats=print_stats)

    if len(subclasses)==0:
        return
    for subclass in subclasses:
        catsubclass = cat[cat[subclass_col] == subclass]
        print(f'# events of {subclass} = {len(catsubclass)}')
        for i,col1 in enumerate(mag_columns):
            for i2, col2 in enumerate(mag_columns):
                if col1==col2:
                    continue
                try:
                    print(' ')
                    #[c, m, r2, fig, ax] = linregress(catsubclass, col1, col2, plot=plot, print_stats=print_stats)
                    if mfixed:
                        fix_slope(catsubclass, col1, col2, mfixed=mfixed[i2]/mfixed[i], plot=plot, print_stats=print_stats)
                    else:
                        linregress(catsubclass, col1, col2, plot=plot, print_stats=print_stats)                    
                except Exception as e:
                    print(f'Failed linregress {col2} vs {col1}')
                    print(e)
                    
def linearRegressMagnitudesBySubclass2ways(cat, subclass_col='subclass', mag_columns=['ML', 'ME'], \
                                      subclasses=['r', 'e', 'l', 'h', 't'], \
                                      plot=False, print_stats=False):



    print(f'# events = {len(cat)}')
    for i,col1 in enumerate(mag_columns):
        for i2,col2 in enumerate(mag_columns[i+1:]): # alt is for col2 in mag_columns[i+1:]:
            if not col1==col2:
                linregress_2ways(cat, col1, col2, plot=False, print_stats=False) 

    if len(subclasses)==0:
        return
    for subclass in subclasses:
        catsubclass = cat[cat[subclass_col] == subclass]
        print(f'# events of {subclass} = {len(catsubclass)}')
        for i,col1 in enumerate(mag_columns):
            for i2, col2 in enumerate(mag_columns):
                if col1==col2:
                    continue
                try:
                    print(' ')
                    try:
                        linregress_2ways(catsubclass, col1, col2, plot=False, print_stats=False)     
                    except:
                        pass
                except Exception as e:
                    print(f'Failed linregress {col2} vs {col1}')
                    print(e)                    


def linregress_2ways(catdf, col1, col2, plot=False, print_stats=False, r2min=.5):
    df = catdf.copy()
    df = df[[col1, col2]]
    df[col1]=df[col1].astype(float)
    df[col2]=df[col2].astype(float)
    df.replace([np.inf, -np.inf, 'NaN'], np.nan, inplace=True)
    df.dropna()
    
    [c1, m1, r21, fig, ax] = linregress(df, col1, col2, plot=False, print_stats=print_stats, show_now=False)
    [c2, m2, r22, fig, ax] = linregress(df, col2, col1, plot=False, print_stats=print_stats, show_now=False) 
    if not m2:
        return
    m = np.mean([m1, 1/m2])
    c = np.mean([c1, -c2/m2])
    r2 = np.mean([r21, r22])
    if r2<r2min:
        return
    print('\nAverage relationship:')
    if c >= 0.0:
        print(f"{col2} = {sigfigs(m,n=3)} {col1} + {sigfigs(c,n=3)}")
        print(f"alt: {col1} = {sigfigs(1/m,n=3)} {col2} - {sigfigs(c/m,n=3)}")       
    else:
        print(f"{col2} = {sigfigs(m,n=3)} {col1} - {sigfigs(c*-1,n=3)}")   
        print(f"alt: {col1} = {sigfigs(1/m,n=3)} {col2} + {sigfigs(c*-1/m,n=3)}")   
    y_pred =  df[col1]*m + c
    print(f"R-squared value: {sigfigs(r2,n=3)}")
    #print("Mean squared error: %.2f" % mean_squared_error(df[col2], y_pred))
    fig, ax = plt.subplots()

    plt.scatter(df[col1], df[col2], color="black")
    plt.plot(df[col1], y_pred, color="blue", linewidth=3)
    plt.xlabel(col1)
    plt.ylabel(col2)
    plt.show()

def fix_slope(df, col1, col2, mfixed=None, cfixed=None, plot=False, print_stats=True):
    [c, m, r2, fig, ax] = linregress(df, col1, col2, plot=plot, print_stats=print_stats, show_now=False)
    if mfixed:
        df[col1]=df[col1].astype(float)
        df[col2]=df[col2].astype(float)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)       
        cf = df[col2] - mfixed * df[col1]
        cfmean = sigfigs(cf.mean(),3)
        
        
        df2 = df.copy()
        col3 = col2+'_fit'
        col4 = col2+'_fixed'
        df2[col3] =  df2[col1]*m + c
        df2[col4] =  df2[col1]*mfixed + cfmean
        if cfmean >= 0.0:
            print(f"{col4} = {sigfigs(mfixed,n=3)} {col1} + {sigfigs(cfmean,n=3)}")
            print(f"alt: {col1} = {sigfigs(1/mfixed,n=3)} {col4} - {sigfigs(cfmean/mfixed,n=3)}")
        else:
            print(f"{col4} = {sigfigs(mfixed,n=3)} {col1} - {sigfigs(cfmean*-1,n=3)}")       
            print(f"alt: {col1} = {sigfigs(1/mfixed,n=3)} {col4} + {sigfigs(cfmean*-1/mfixed,n=3)}") 
        
        fig, ax = plt.subplots()
        #xmin=df2[col1].min()
        #xmax=df2[col2].max()
        #x=np.arange(xmin,xmax,(xmax-xmin)/100)
        #y1 = x*m + c
        #y2 = x*mfixed + cf.mean()
        df2.plot.scatter(x=col1, y=col2, rot=90, ax=ax)
        df2.plot.line(x=col1, y=col3, rot=90, ax=ax)
        df2.plot.line(x=col1, y=col4, rot=90, ax=ax)  
        #plt.plot(x,y1)
        #plt.plot(x,y2)        
        #df2.plot.line(x=col1, y=col4, rot=90, ax=ax) 
        plt.show()

def local_magnitude(A, R, n=1.00, k=3.01e-6, c=0.70, roundDigits=2):
    '''
    ML = log A + n log R + k R + c

    Bakun & Joyner [1984] found n = 1.0 (1/R decay), k = 3.01e-6, c = 0.70 (after I convereted to S.I. units)
    Hutton & Boore [198?] found n = 1.11, k = 1.89e-6, c = 0.261

    A and R must be in m (in papers, A in mm, R in km, but I converted to S.I. units)

    ML defined as ML=3 for A=1e-3 m at 1e5 m (1 mm at 100 km)
    Hutton & Boore redefine as ML=3 for A=1e-2 m at 1.7e4 m (10 mm at 17 km)
    Both scales give these results.

    '''
    ML = np.log10(A) + n * np.log10(R) + k * R + c
    if roundDigits:
        return np.round(ML, roundDigits)
    else:
        return ML

def seismic_energy(eng, R, rho_e=2500, c_e=2500, fs=75):
    E = 2 * np.pi * R**2 * rho_e * c_e * eng / fs
    return E

def energy_magnitude(E, a=-3.2, b=2/3, roundDigits=2):
    # should a=-3.2?
    ME = b * np.log10(E) + a
    if roundDigits:
        return np.round(ME, roundDigits)
    else:
        return ME   
        
##### functions from 510_Instrumental_VEI

# Load raw seismic data - and set units accordingly
DATA_DIR = os.path.join('data')
SDS_DIR = os.path.join(DATA_DIR, 'continuous','SDS')
SAM_DIR = os.path.join(DATA_DIR, 'continuous','SAM')
RESPONSE_DIR = os.path.join(DATA_DIR, 'responses')

# Load raw seismic data - and set units accordingly
SDS_DIR = os.path.join('/data', 'SDS')
#SAM_DIR = os.path.join(DATA_DIR, 'continuous','SAM')
RESPONSE_DIR = os.path.join('data', 'responses')

#resultsDF = pd.DataFrame(columns=['Event', 'sum(ER)', 'ME', 'sum(ER_VLP)', 'ME_VLP', 'DR', 'DR_VLP', 'DRS', 'DRS_VLP', 'M_DR', 'M_DR_VLP', 'M_DRS', 'M_DRS_VLP'])

resultsDF = pd.DataFrame(columns=['Event', 'start', 'end', 'duration', 'ML', 'sum(ER)', 'ME', 'DR', 'DRS',  'M_DR',  'M_DRS'])

def magnitude_DR(v, return_dataframe=True, correction=0.1):

    '''
    Based on:
        ML = log10(A) + 1.11 log10(R) -1.89e-6 R + 3.58
    which is Hutton & Boore [1987] in SI units
    
    Assuming a sine wave:
        DR = rms(A).R = A/sqrt(2) . R = A/1.414 . R
    Hence:
        ML = log10(1.414 DR/R) + 1.11 log10(R) - 1.89e-6 R + 3.58
        ML = 0.15 + log10(DR) + 0.11 log10(R) - 1.89e-6 R + 3.58
        ML = log10(DR) + 0.11 log10(R) - 1.89e-6 R + 3.73
    For R = 1000 m, 
        ML = log10(DR) + 4.06
    For R = 10,000 m:
        ML = log10(DR) + 4.15
    Therefore, for most volcano-seismic stations, 1-10 km:
        ML = log10(DR) + 4.1
    However, DR here is expressed in m^2. If we convert to cm^2:
        ML = log10(DR) + 0.1
    '''

    if isinstance(v, float):
        m = np.log10(v) + correction
        return round(m, 1)
    elif isinstance(v, DR) or isinstance(v, DRS):
        thisDRobj = v
        ids = thisDRobj.get_seed_ids()
    
        lod = []
        lom = []
        
        for id in ids:
            maxDR = thisDRobj.dataframes[id]['rms'].max() 
            M_DR = np.log10(maxDR) + correction
            #print(f'id={id}, maxDR={maxDR:.1f}, M_DR={M_DR:.1f}')
            lod.append({'id':id, 'maxDR':round(maxDR,1), 'M_DR':round(M_DR,1) }) 
            lom.append(M_DR)
        
        if return_dataframe:
            df = pd.DataFrame.from_dict(lod)
            return df
        else:
            return round(np.nanmedian(np.array(lom)),1)
    else:
        return None


def compute_ML(stD, inv, source, f=1.0):
    # from: https://gfzpublic.gfz-potsdam.de/rest/items/item_43489/component/file_56102/content
    # find max of each trace
    ML = []
    for tr in stD:
        tr2 = tr.copy()
        tr2.filter('bandpass', freqmin=0.5, freqmax=10.0, corners=2) # similar frequency range to Wood-Anderson

        # compute R
        R = get_distance_km(tr2.id, inv, source) * 1000
        
        #A = np.nanmax(tr2.data) * 1e9 # amplitude in nm
        #ML.append(np.log10(A) + 1.11 * np.log10(R) + 0.00189 * R - 2.09)

        # SI units equivalent
        A = np.nanmax(tr2.data) # Amplitude in m
        #thisML = np.log10(A) + f * np.log10(R) + 1.89e-6*R + 3.58, f=1.11
        
        #thisML = np.log10(A) + f * np.log10(R) + 3.01e-6*R + 0.70
        thisML = local_magnitude(A, R)
        print(f'A={A}, R={R}, ML={thisML}')
        ML.append(thisML)
    
    return round(np.median(np.array(ML)), 1)
    

from obspy.geodetics.base import gps2dist_azimuth, degrees2kilometers
def get_distance_km(seed_id, inventory, source):
    coordinates = {}
    coordinates[seed_id] = inventory.get_coordinates(seed_id)
    if seed_id[0:2]=='MV':
        if coordinates[seed_id]['longitude'] > 0:
            coordinates[seed_id]['longitude']  *= -1
    distance_m, az_source2station, az_station2source = gps2dist_azimuth(source['lat'], source['lon'], coordinates[seed_id]['latitude'], coordinates[seed_id]['longitude'])
    distance_km = distance_m/1000
    return distance_km

# Load data from SDS, convert into raw stream object
from obspy.clients.filesystem.sds import Client as sdsclient
#from obspy.signal.freqattributes import welch

from scipy.signal import welch
def spectralRatio(tr, paddingTime):

    print(tr)

    nfft = 8192
    fs = tr.stats.sampling_rate
    df = fs/nfft
    f = np.arange(0.0, fs, df)
    
    noise = tr.copy().trim(endtime=tr.stats.starttime+paddingTime)
    f, nPxx_spec = welch(noise.data, fs, nperseg=nfft, scaling='spectrum')
    '''
    plt.figure()
    plt.semilogx(f, np.sqrt(nPxx_spec))
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Linear NOISE spectrum [RMS]')
    plt.xlim([0.01, 27.0])
    plt.show()
    '''

    signal = tr.copy().trim(starttime=tr.stats.starttime+paddingTime, endtime=tr.stats.endtime-paddingTime )
    f, sPxx_spec = welch(signal.data, fs, nperseg=nfft, scaling='spectrum')
    '''
    plt.figure()
    plt.semilogx(f, np.sqrt(sPxx_spec))
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Linear SIGNAL spectrum [RMS]')
    plt.xlim([0.01, 27.0])
    plt.show()
    '''

    plt.figure()
    plt.semilogx(f, np.sqrt(sPxx_spec)/np.sqrt(nPxx_spec))
    plt.xlim([0.01, 27.0])
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Linear SIGNAL/NOISE spectrum [RMS]')
    plt.show()  

    print('\n\n\n')
 

import matplotlib.pyplot as plt
def sds2Stream(SDS_DIR, startTime, endTime, network='MV', station='*', location='*', channel="[SBEHCD]H*", padding=0.2, \
               inv=None, pre_filt = [0.02, 0.03, 18, 27], miniseed=True, trace_ids=None, spectralRatioPlots=False ):
    mySDSclient = sdsclient(SDS_DIR)
    paddingTime = (endTime-startTime)*padding
    if paddingTime < 120:
        paddingTime = 120

    print(f'Loading data from SDS archive for {network}.{station}.{location}.{channel} from {startTime} to {endTime} with {paddingTime} padding seconds')
    if trace_ids:
        stR = obspy.Stream()
        for id in trace_ids:
            net, sta, loc, cha = id.split('.')
            this_st = mySDSclient.get_waveforms(net, sta, loc, cha, startTime-paddingTime, endTime+paddingTime)
            for tr in this_st:
                stR.append(tr)
    else:
        stR=mySDSclient.get_waveforms(network, station, location, channel, startTime-paddingTime, endTime+paddingTime)
    for tr in stR:
        if np.all(tr.data==tr.data[0]):
            # all values same
            print(f'Flat waveform: removing {tr.id}')
            stR.remove(tr)

    if len(stR)==0:
        print('No data - will try loading miniseed files directly')
        print(mySDSclient)
        nslc_list = mySDSclient.get_all_nslc(datetime=startTime)
        for nslc in nslc_list:
            filename = mySDSclient._get_filename(nslc[0], nslc[1], nslc[2], nslc[3], startTime)
            if os.path.isfile(filename):
                yes = 1
            else:
                yes = 0
            pct, ngaps = mySDSclient.get_availability_percentage(nslc[0], nslc[1], nslc[2], nslc[3], startTime, endTime)
            print(f"{nslc}, {filename}, exists={yes}, availability={pct}")
            if yes:
                os.system(f'ls -l {filename}')
                thisst = obspy.read(filename, 'MSEED')
                thisst.merge()
                for tr in thisst:
                    stR.append(tr)
        print('Full day in SDS archive:')
        print(stR)
        stR.trim(starttime=startTime, endtime=endTime)
        for tr in stR:
            if tr.stats.npts==0:
                stR.remove(tr)
        if len(stR)==0:
            del mySDSclient
            raise Exception("no data loaded for this event")
            return

    del mySDSclient    
    print('raw seismograms')
    print(stR)
    #stR.plot(equal_scale=False);

    if spectralRatioPlots:
        for tr in stR:
            try:
                spectralRatio(tr, paddingTime)
            except:
                pass

    if inv:
        print('\nCorrecting to velocity')
        stV = stR.copy().detrend('linear')
        for tr in stV: # could change channel filter above and add code here to apply overall sensitivity correction to infrasound
            try:
                if tr.stats.channel[0] in 'ES':   
                    tr.remove_response(inventory=inv, pre_filt=[0.33, 0.5, pre_filt[2], pre_filt[3]], output="VEL", plot=False, water_level=None)  
                else:
                    tr.remove_response(inventory=inv, pre_filt=pre_filt, output="VEL", plot=False, water_level=None)                      
                print(f'{tr.id} corrected')
            except Exception as e:
                print(e)
                print(f'removing {tr.id} from stV')
                stV.remove(tr)
        stV.trim(starttime=startTime, endtime=endTime)
        print('velocity seismograms')
        print(stV)
        #stV.plot(equal_scale=False);

        print('\nCorrecting to displacement')
        stD = stR.copy().detrend('linear')
        for tr in stD:
            try:
                if tr.stats.channel[0] in 'ES': 
                    tr.remove_response(inventory=inv, pre_filt=[0.33, 0.5, pre_filt[2], pre_filt[3]], output="DISP", plot=False, water_level=None)  
                else:
                    tr.remove_response(inventory=inv, pre_filt=pre_filt, output="DISP", plot=False, water_level=None)   
                print(f'{tr.id} corrected')
            except Exception as e:
                print(e)
                print(f'removing {tr.id} from stD')
                stD.remove(tr)
        stD.trim(starttime=startTime, endtime=endTime)
        print('displacement seismograms')
        print(stD)
        #stD.plot(equal_scale=False);

        if miniseed:
            if not os.path.isdir('data4Felix'):
                os.makedirs('data4Felix')
            mseedfile = os.path.join('data4Felix', f"{network}_{startTime}_{endTime}.mseed")
            print(f'Saving displacement seismograms to {mseedfile}')
            stD.write(mseedfile, format='mseed')
            inv.write(mseedfile.replace('.mseed', '.xml'), format='stationxml')
        return stV, stD
    else:
        return stR


def compute_sam_metrics(stV, stD, sampling_interval=60, filter=[0.5, 18.0], \
                        bands={'VLP': [0.1, 0.5], 'LP':[0.5, 4.0], 'VT':[4.0, 18.0]}, corners=2):
    print('computing DSAM and VSEM')
    dsamObj = DSAM(stream=stD, sampling_interval=sampling_interval, filter=filter, bands=bands, corners=corners)
    vsemObj = VSEM(stream=stV, sampling_interval=sampling_interval, filter=filter, bands=bands, corners=corners)
    return vsemObj, dsamObj

def compute_reduced_metrics(eventname, stime, etime, dsamObj, vsemObj, source, inv, resultsDF, ML, surfaceWaveSpeed_kms=2.0, peakf=None):

    dsamObj = dsamObj.select(component='Z')

    # reduced energy
    ERobj = vsemObj.compute_reduced_energy(inv, source)
    #ERobj.plot(metrics='energy')
    sumE, ME = ERobj.sum_energy(metric='energy') # for Choy & Boatright b=2/3 a=-3.2

    '''
    ERobjBB = ERobj.select(channel='[BH]HZ')
    if len(ERobjBB)>0:
        sumE_VLP, ME_VLP = ERobjBB.sum_energy(metric='VLP')
    else:
        sumE_VLP = 0
        ME_VLP = -np.Inf
    '''
    
    # body wave reduced displacement
    print('Computing body-wave DR')
    DRobj = dsamObj.compute_reduced_displacement(inv, source, surfaceWaves=False, Q=None, wavespeed_kms=surfaceWaveSpeed_kms, peakf=peakf)
    #DRobj.plot(metrics='rms')
    DRmaxrms = sigfigs(DRobj.max(metric='rms'), n=2)
    #M_DR = magnitude_DR(DRmaxrms)

    '''
    print('Computing body-wave DR in the VLP band')
    DRobjBB = DRobj.select(channel='[BH]HZ')
    print(DRobj.get_seed_ids(), DRobjBB.get_seed_ids() )1.89*1.7
    if len(DRobjBB)>0:
        DRmaxvlp = sigfigs(DRobjBB.max(metric='VLP'), n=2)
        M_DR_vlp = magnitude_DR(DRmaxvlp)
    else:
        DRmaxvlp  = 0
        M_DR_vlp = -np.Inf 
    '''
        
    # surface wave reduced displacement
    print('Computing surface-wave DRS')
    DRSobj = dsamObj.compute_reduced_displacement(inv, source, surfaceWaves=True, Q=None, wavespeed_kms=surfaceWaveSpeed_kms, peakf=peakf)
    #DRSobj.plot(metrics='rms')   
    DRSmaxrms = sigfigs(DRSobj.max(metric='rms'), n=2)
    #M_DRS = magnitude_DR(DRSmaxrms)
    #M_DRS = magnitude_DR(DRSobj, return_dataframe=False)

    '''
    print('Computing surface-wave DRS in the VLP band')
    DRSobjBB = DRSobj.select(channel='[BH]HZ')
    if len(DRSobjBB)>0:
        DRSmaxvlp = sigfigs(DRSobjBB.max(metric='VLP'), n=2)
        M_DRS_vlp = magnitude_DR(DRSmaxvlp)
        #M_DRS_vlp = magnitude_DR(DRSobjBB, return_dataframe=False)
    else:
        #sumE_VLP, ME_VLP = ERobj.sum_energy(metric='VLP')
        DRSmaxvlp  = 0
        M_DRS_vlp = -np.Inf
    '''

    # add to results
    '''
    resultsDF.loc[len(resultsDF.index)] = [eventname, sigfigs(sumE,n=2), ME, sigfigs(sumE_VLP,n=2), ME_VLP, DRmaxrms, DRmaxvlp, \
                                           DRSmaxrms, DRSmaxvlp, M_DR, M_DR_vlp, M_DRS, M_DRS_vlp]
    '''
    duration = etime-stime
    resultsDF.loc[len(resultsDF.index)] = [eventname, stime.strftime('%Y-%m-%d %H:%M'), etime.strftime('%Y-%m-%d %H:%M'), \
                                           round(duration,0), ML, sigfigs(sumE,n=2), round(ME,1), DRmaxrms, DRSmaxrms] #, M_DR, M_DRS]
        

def wrapper(SDS_DIR, network, startTime, endTime, stationxml, resultsDF, eventname, source, sampling_interval=60, \
            trace_ids=None, pre_filt=[0.4, 0.6, 13, 20], filter=[0.5, 18.0], \
                        bands={'VLP': [0.1, 0.5], 'LP':[0.5, 4.0], 'VT':[4.0, 18.0]}, corners=2):

    inv = obspy.read_inventory(stationxml)

    if isinstance(startTime, list):
        for i in range(len(startTime)):
            if isinstance(eventname, list):
                subeventname = eventname[i]
            else:
                subeventname = f"{eventname}: subevent {i}"
            wrapper(SDS_DIR, network, startTime[i], endTime[i], stationxml, resultsDF, subeventname, source, sampling_interval=sampling_interval, trace_ids=trace_ids)

    try:
        stV, stD = sds2Stream(SDS_DIR, startTime, endTime, network=network, station='*', location='*', channel="[SBEHCD]H*", padding=0.1, \
               inv=inv, pre_filt = pre_filt, miniseed=True, trace_ids=trace_ids )
    except Exception as e:
        print(e)
        #print(type(startTime), type(endTime))
        #print(inv)
        return
    else:
        if len(stV)>0:
            ML = compute_ML(stD, inv, source, f=1.0)
        
            delta = sampling_interval
            if endTime-startTime < 7200:
                delta = (endTime-startTime)/120
            
            vsemObj, dsamObj = compute_sam_metrics(stV, stD, sampling_interval=delta, filter=filter, \
                                bands=bands, corners=corners)
        
            compute_reduced_metrics(eventname, startTime, endTime, dsamObj, vsemObj, source, inv, resultsDF, ML, surfaceWaveSpeed_kms=2.0, peakf=None)
        
            #display(resultsDF)    