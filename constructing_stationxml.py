def load_mvo_inventory(tr, CALDIR):
    this_inv = None
    matchcode = None
    if len(tr.stats.channel)<3:
        return this_inv
    if tr.stats.channel[0] in 'ES':
        matchcode = '[ES]'
    elif tr.stats.channel[0] in 'BH':
        matchcode = '[BH]'
    if not matchcode:
        print("Cannot match trace ID %s ",tr.id)
        return this_inv
    #xmlfilepattern = os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) )
    xmlfiles = glob(os.path.join(CALDIR, "station.MV.%s..%s*%s.xml" % (tr.stats.station, matchcode, tr.stats.channel[2]) ))
    N = len(xmlfiles)
    if N==1:
        xmlfile = xmlfiles[0]
        print('Correcting %s with %s' % (tr.id, xmlfile))
        this_inv = read_inventory(xmlfile)  
    return this_inv

def respfiles2masterstationxml(SEISAN_DATA, xmlfile):# Merge Montserrat RESP files
    # A function we are only likely to use once, but an important one to keep
    from obspy.io.xseed.core import _read_resp
    from lib.libseisGT_old3 import create_dummy_inventory, merge_inventories
    station0hypfile = os.path.join(SEISAN_DATA, 'DAT', 'STATION0_MVO.HYP')
    station_locationsDF = parse_STATION0HYP(station0hypfile) 
    respfiles = glob.glob(os.path.join(SEISAN_DATA, 'CAL', 'RESP*'))
    master_inv = create_dummy_inventory()
    for respfile in respfiles:
        this_inv = None
        print('RESP file = ',respfile)
        this_inv = _read_resp(respfile)
        this_inv = inventory_fix_id_mvo(this_inv)
        sta_code = this_inv.networks[0].stations[0].code
        location = station_locationsDF[station_locationsDF['name']==sta_code].iloc[0].to_dict()
        for station in this_inv.networks[0].stations:
            station.latitude = location['lat']
            station.longitude = location['lon']
            station.elevation = location['elev']
            for channel in station.channels:
                channel.latitude = location['lat']
                channel.longitude = location['lon']
                channel.elevation = location['elev']
                channel.depth = 0.0
                if channel.sample_rate==75.19:
                    channel.sample_rate=75.0
        merge_inventories(master_inv, this_inv)  
    master_inv.write(xmlfile,format="STATIONXML", validate=True)