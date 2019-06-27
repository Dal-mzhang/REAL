import obspy
import os
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
build_taup_model("../tt_db/itvel.nd")
model = TauPyModel(model="itvel")

if not os.path.exists('pk'):
    os.makedirs('pk')

with open("../station.dat", "r") as f:
    for station in f:
        stlo, stla, net, sta, chan, elev = station.split()
        print(net,sta)
        output1 = './pk/'+net+'.'+sta+'.'+'P.txt'
        f1 = open(output1,'w')
        output2 = './pk/'+net+'.'+sta+'.'+'S.txt'
        f2 = open(output2,'w')

        with open("catalog_use.txt","r") as p:
            for event in p:
                time,lon,lat,dep,mag = event.split()
                date,time = time.split('T')
                hour,minute,sec0 = time.split(':')
                sec = sec0[:5]
                torg = float(hour)*3600 + float(minute)*60 + float(sec)
                ddist,az1,az2 = gps2dist_azimuth(lat1=float(lat),lon1=float(lon),lat2=float(stla),lon2=float(stlo))
                ddist = (ddist/1000)/111.19
                #print(ddist,dep)
                arrivals = model.get_travel_times(source_depth_in_km=float(dep), distance_in_degree=ddist, phase_list=["P","p","S","s"])
                i = 0
                pi = 0
                si = 0
                while(i<len(arrivals)):
                    arr = arrivals[i]
                    i = i + 1
                    if((arr.name == 'P' or arr.name == 'p') and pi == 0):
                        pname = arr.name
                        tp = torg + arr.time
                        pi = 1
                        f1.write('{} 1.0 0.0\n'.format(tp))

                    if((arr.name == 'S' or arr.name == 's') and si == 0):
                        sname = arr.name
                        ts = torg + arr.time
                        si = 1
                        f2.write('{} 1.0 0.0\n'.format(ts))
                    
                    if(pi == 1 and si == 1):
                        break
        f1.close()
        f2.close()
