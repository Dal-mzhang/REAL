from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import os

source="INGV"
stationfile="../ingv.sta"

year = 2016
mon = 10
day = 14

tsec = 86400

dirdata = str(year)+str(mon)+str(day)+'/'
tb = str(year)+'-'+str(mon)+'-'+str(day)+"T00:00:00.00Z"

if not os.path.exists(dirdata):
    os.makedirs(dirdata)

client = Client(source)
tbeg=UTCDateTime(tb)
tend = tbeg + tsec

with open(stationfile, "r") as f:
    for station in f:
        stlo, stla, net, sta, chan, elev = station.split()
        chane = chan[:2]+"E"
        chann = chan[:2]+"N"
        chanz = chan[:2]+"Z"
        chan0 = [chane,chann,chanz]
        for chan1 in chan0:
            print(sta,chan1)
            try:
                st = client.get_waveforms(network=net,station=sta,channel=chan1,starttime=tbeg,endtime=tend,location=False,attach_response=True)
                st.merge(method=1, fill_value='interpolate')
                st.interpolate(sampling_rate=100,startime=tbeg) # like sac interpolate
                st.detrend("demean")
                st.detrend("linear")
                pre_filt = [0.001, 0.002, 25, 30] # bandpass filter
                st.remove_response(pre_filt=pre_filt,water_level=10,taper=True,taper_fraction=0.00001)
                st.write(filename=dirdata+net+'.'+sta+'.'+chan1+".SAC",format="SAC")
            except:
                print("doesn't exist",sta,chan1)
