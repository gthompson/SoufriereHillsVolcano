from obspy.clients.nrl import NRL
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy import UTCDateTime, read_inventory
from obspy.core.inventory.response import Response, InstrumentSensitivity #PolesZerosResponseStage

from libseisGT import add_to_trace_history


import os
import requests
from obspy.io.xseed import Parser

""" the idea of this is to give ways to just use the overall sensitivyt easily, or the full
 instrument response for different types of datalogger/sensors USF has used """

def get_rboom():
    return get_rsb()

def get_rsb():
    # Step 1: Download the StationXML file from the provided URL
    url = 'https://manual.raspberryshake.org/_downloads/57ab6152abedf7bb15f86fdefa71978c/RSnBV3-20s.dataless.xml-reformatted.xml'
    xmlfile = 'rsb_v3.xml'
    if not os.path.isfile(xmlfile):
        response = requests.get(url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Save the content to a local file
            with open(xmlfile, "wb") as f:
                f.write(response.content)
            print("File downloaded successfully.")

        else:
            print(f"Failed to download file. HTTP status code: {response.status_code}")
            return None

    # Step 2: Read the downloaded StationXML file into an ObsPy Inventory object
    try:
        inventory = read_inventory(xmlfile)
    except:
        responses = {}
        poles = [-0.312 + 0.0j, -0.312 + 0.0j]
        zeros = [0.0 + 0j, 0.0 + 0j]
        sensitivity = 56000
        # note: cannot use Pa in response, so use m/s but remember it is really Pa
        responses['HDF'] = Response.from_paz(zeros, poles, sensitivity, input_units='m/s', \
                                                         output_units='Counts' )
        poles = [-1.0, -3.03, -3.03 ]
        zeros = [0.0 + 0j, 0.0 + 0j, 0.0 + 0j]
        sensitivity = 399650000
        responses['EHZ'] = Response.from_paz(zeros, poles, sensitivity, input_units='m/s', \
                                                         output_units='Counts' )
        return responses

    return inventory

def get_rs1d_v4():
    # Step 1: Download the StationXML file from the provided URL
    url = 'https://manual.raspberryshake.org/_downloads/e324d5afda5534b3266cd8abdd349199/out4.response.restored-plus-decimation.dataless'
    
    responsefile = 'rs1d_v4.dataless'
    xmlfile = 'rs1d_v4.xml'
    if not os.path.isfile(responsefile):
        response = requests.get(url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Save the content to a local file
            with open(responsefile, "wb") as f:
                f.write(response.content)
            print("File downloaded successfully.")

        else:
            print(f"Failed to download file. HTTP status code: {response.status_code}")
            return None
    sp = Parser(responsefile)
    sp.write_seed(xmlfile)

    # Step 2: Read the downloaded StationXML file into an ObsPy Inventory object
    inventory = read_inventory(xmlfile)
    return inventory

def get_rs3d_v5():
    # Step 1: Download the StationXML file from the provided URL
    url = 'https://manual.raspberryshake.org/_downloads/858bf482b1cd1b06780e9722c1d0b2db/out4.response.restored-EHZ-plus-decimation.dataless-new'
    
    responsefile = 'rs3d_v5.dataless'
    xmlfile = 'rs3d_v5.xml'
    if not os.path.isfile(responsefile):
        response = requests.get(url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Save the content to a local file
            with open(responsefile, "wb") as f:
                f.write(response.content)
            print("File downloaded successfully.")

        else:
            print(f"Failed to download file. HTTP status code: {response.status_code}")
            return None
    sp = Parser(responsefile)
    sp.write_seed(xmlfile)

    # Step 2: Read the downloaded StationXML file into an ObsPy Inventory object
    inventory = read_inventory(xmlfile)

    return inventory

def get_rs3d_v3():
    # Step 1: Download the StationXML file from the provided URL
    url = 'https://manual.raspberryshake.org/_downloads/85ad0fd97072077ea8a09e42ab15db4b/out4.response.restored-plus-decimation.dataless'
    
    responsefile = 'rs3d_v3.dataless'
    xmlfile = 'rs3d_v3.xml'
    if not os.path.isfile(responsefile):
        response = requests.get(url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Save the content to a local file
            with open(responsefile, "wb") as f:
                f.write(response.content)
            print("File downloaded successfully.")

        else:
            print(f"Failed to download file. HTTP status code: {response.status_code}")
            return None
    sp = Parser(responsefile)
    sp.write_seed(xmlfile)

    # Step 2: Read the downloaded StationXML file into an ObsPy Inventory object
    inventory = read_inventory(xmlfile)
    return inventory

def centaur(inputVoltageRange):
    countsPerVolt = 0.4e6 * 40/inputVoltageRange;
    return countsPerVolt

def trillium():
    voltsPerMS = 754; # V / (m/s)
    return voltsPerMS

def infraBSU(HgInThisSensor=0.5):
    # model 0.5" is default
    # 0.1 - 40 Hz flat
    oneInchHg2Pa = 3386.4;
    linearRangeInPa = oneInchHg2Pa * HgInThisSensor;
    selfNoisePa = 5.47e-3;
    voltsPerPa = 46e-6; # from infraBSU quick start guide
    return voltsPerPa

def ChaparralM25():
    # 36 V p2p
    # 0.1 - 200 Hz flat
    selfNoisePa = 3e-3;
    voltsPerPaHighGain = 2.0; # 18 Pa full scale. 
    voltsPerPaLowGain = 0.4; # 90 Pa full scale. recommended for 24-bit digitizers. 
    voltsPerPaVMod = 0.4 * 90/720; # 720 Pa full scale.
    # Volcano mod reduces sensitivity further.
    return voltsPerPaLowGain

countsPerMS = centaur(40.0) * trillium()
countsPerPa40 = centaur(40.0) * infraBSU(0.5)
countsPerPa1 = centaur(1.0) * infraBSU(0.5)
countsPerPaChap40 = centaur(40.0) * ChaparralM25()
countsPerPaChap1 = centaur(1.0) * ChaparralM25()
rs1d    = 469087000 # counts/m/s
rboom   = 56000 # counts/Pa
rsb     = 399650000 # counts/m/s for a shake&boom EHZ
rs3d    = 360000000

def correctUSFstations(st, apply_calib=False, attach=True):
    chap25 = {}
    datalogger = 'Centaur'
    sensor = 'TCP'
    Vpp = 40
    inventories = Inventory()
    for tr in st:
        if tr.stats.network=='AM':
            continue
        if tr.stats.sampling_rate>=50:
            #tr.detrend()
            if tr.stats.channel[1]=='H': # a seismic velocity high-gain channel. L for low gain, N for accelerometer
                # trace is in counts. there are 0.3 counts/ (nm/s).
                #tr.data = tr.data / 0.3 * 1e-9 # now in m/s
                calib = countsPerMS
                units = 'm/s'
                                
            if tr.stats.channel[1]=='D': # infraBSU channel?
                #trtr.data = tr.data / countsPerPa40
                # Assumes 0.5" infraBSU sensors at 40V p2p FS
                # But could be 1" or 5", could be 1V p2p FS or could be Chaparral M-25
                units = 'Pa'
                sensor = 'infraBSU'
                Vpp = 1
                calib = countsPerPa1 # default is to assume a 0.5" infraBSU and 1V p2p
                # First let's sort out all the SEED ids that correspond to the Chaparral M-25
                if tr.stats.station=='BCHH1' and tr.stats.channel[2]=='1':
                    calib = countsPerPaChap40
                    chap25 = {'id':tr.id, 'p2p':40}
                    sensor = 'Chaparral'
                    Vpp = 40
                elif tr.id == 'FL.BCHH3.10.HDF':
                    if tr.stats.starttime < UTCDateTime(2022,5,26): # Chaparral M25. I had it set to 1 V FS. Should have used 40 V FS. 
                        calib = countsPerPaChap1
                        chap25 = {'id':tr.id, 'p2p':1}
                        sensor = 'Chaparral'
                        Vpp = 1
                    else:
                        calib = countsPerPaChap40
                        chap25 = {'id':tr.id, 'p2p':40}
                        sensor = 'Chaparral'
                        Vpp = 40
                elif tr.id=='FL.BCHH4.00.HDF':
                    calib=countsPerPaChap40    
                    sensor = 'Chaparral'
                    Vpp = 40
                # anything left is infraBSU and we assume we used a 0.5" sensor, but might sometimes have used the 1" or even the 5"
                # assume we used a 1V peak2peak, except for the following cases
                elif tr.stats.station=='BCHH' or tr.stats.station=='BCHH1' or tr.stats.station[0:3]=='SAK' or tr.stats.network=='NU':
                    calib = countsPerPa40
                    Vpp = 40
                elif chap25:
                        net,sta,loc,chan = chap25['id'].split()
                        if tr.stats.network == net and tr.stats.station == sta and tr.stats.location == loc:
                            if chap25['p2p']==40:
                                calib = countsPerPa40
                                Vpp = 40
                            else:
                                Vpp = 1
            inv = make_inv(tr.stats.network, tr.stats.station, tr.stats.location, tr.stats.channel, \
                           datalogger=datalogger, sensor=sensor, Vpp=Vpp, fsamp=tr.stats.sampling_rate, \
                              sensitivty=calib, units=units)
            inventories = inventories + inv
            if apply_calib:
                tr.data = tr.data / calib
                tr.stats['sensitivity'] = calib
                tr.stats['units'] = units
                add_to_trace_history(tr, 'calibration_applied')   
            tr.stats['datalogger'] = datalogger
            tr.stats['sensor'] = sensor           
    if attach:
        st.attach_response(inventories)
    return inventories    
     

def make_inv(net, sta, loc, chans, datalogger='Centaur', sensor='TCP', Vpp=40, fsamp=100.0, \
             lat=0.0, lon=0.0, elev=0.0, depth=0.0, sitename='', \
                ondate=UTCDateTime(1970,1,1), offdate=UTCDateTime(2025,12,31), sensitivity=None, units=None):
    nrl = NRL('http://ds.iris.edu/NRL/')
    if datalogger == 'Centaur':
        if Vpp==40:
            datalogger_keys = ['Nanometrics', 'Centaur', '40 Vpp (1)', 'Off', 'Linear phase', "%d" % fsamp]
        elif Vpp==1:
            datalogger_keys = ['Nanometrics', 'Centaur', '1 Vpp (40)', 'Off', 'Linear phase', "%d" % fsamp]
    elif datalogger == 'RT130':
        datalogger_keys = ['REF TEK', 'RT 130 & 130-SMA', '1', "%d" % fsamp]
    else:
        print(datalogger, ' not recognized')
        print(nrl.dataloggers[datalogger])
    print(datalogger_keys)

    if sensor == 'TCP':
        sensor_keys = ['Nanometrics', 'Trillium Compact 120 (Vault, Posthole, OBS)', '754 V/m/s']
    elif sensor == 'L-22':
        sensor_keys = ['Sercel/Mark Products','L-22D','2200 Ohms','10854 Ohms']
    elif sensor[:4] == 'Chap':
        sensor_keys = ['Chaparral Physics', '25', 'Low: 0.4 V/Pa']
    elif sensor.lower() == 'infrabsu':
        #sensor_keys = ['JeffreyBJohnson', 'sensor_JeffreyBJohnson_infraBSU_LP21_SG0.000046_STairPressure', '0.000046 V/Pa']
        sensor_keys = ['JeffreyBJohnson', 'infraBSU', '0.000046 V/Pa']
    else:
        print(sensor, ' not recognized')
        print(nrl.sensors[sensor])
    print(sensor_keys)


    try:
        response = nrl.get_response(sensor_keys=sensor_keys, datalogger_keys=datalogger_keys)
    except:
        #response = nrl.get_response(datalogger_keys=datalogger_keys)
        response = None
        # use just a sensitivity instead

    channels = []
    if not isinstance(chans, list):
        chans = [chans]
    for chan in chans:
        channel = Channel(code=chan,
                      location_code=loc,
                      latitude=lat,
                      longitude=lon,
                      elevation=elev,
                      depth=depth,
                      sample_rate=fsamp,
                      start_date=ondate,
                      end_date=offdate,
                      )
        if response:
            channel.response = response
        elif sensitivity:
            poles=None
            zeros=None
            if not hasattr(channel, "sensitivity"):
                
                if sensor.lower()=='infrabsu':
                    poles = [-0.301593 + 0.0j]
                    zeros = [0.0 + 0j]
                    print(f'Warning: units changed from {units} to m/s because Obspy cannot work with Pa')
                if poles or zeros:
                    # although input units are Pa, obspy can only compute overall sensitivity and remove response if units are seismic.
                    channel.response = Response.from_paz(zeros, poles, sensitivity, input_units='m/s', \
                                                         output_units='Counts' )
                    '''channel.response = PolesZerosResponseStage(1, sensitivity, 1.0, 'Counts', units, pz_transfer_function_type='LAPLACE (RADIANS/SECOND)', \
                                                               normalization_frequency=1.0,  poles=poles, zeros=zeros)
                    '''
                else:
                    sensitivityObj = InstrumentSensitivity(sensitivity, 1.0, units, 'Counts')
                    channel.response = Response(instrument_sensitivity=sensitivityObj)


        channels.append(channel)
    station = Station(code=sta,
                      latitude=lat,
                      longitude=lon,
                      elevation=elev,
                      creation_date=ondate,
                      site=Site(name=sitename),
                      channels=channels,
                      start_date=ondate,
                      end_date=offdate,
                      )

    network = Network(code=net,
                     stations=[station])
    inventory = Inventory(networks=[network], source="USF_instrument_responses.py")
    return inventory




if __name__ == '__main__':
    inv = make_inv('FL', 'BCHH', '', ['HHZ', 'HHN', 'HHE'], datalogger='Centaur', sensor='TCP', Vpp=40, fsamp=100, sitename='Beach House original', ondate=UTCDateTime(2016,2,24) )
    #inv.write("BCHH_original_seismic.sml", format="stationxml", validate=True)
    inv = make_inv('FL', 'BCHH', '', ['HDF'], datalogger='Centaur', sensor='Chap', Vpp=40, fsamp=100, sitename='Beach House Sonic', ondate=UTCDateTime(2017,8,1) )
    #inv.write("BCHH_sonic_Chap.sml", format="stationxml", validate=True) # Replace PA with Pa and 400,000 in InstrumentSensivitity Value with 160,000
   
    print('calibrations:')
    print('trillium+centaur40 = %f' % countsPerMS)        
    print('infraBSU+centaur40 = %f' % countsPerPa40)        
    print('infraBSU+centaur1 = %f' % countsPerPa1)        
    print('chaparralM25+centaur40 = %f' % countsPerPaChap)        
    print('chaparralM25+centaur1 = %f' % (countsPerPaChap*40))        
