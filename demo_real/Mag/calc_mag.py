#calculate local magnitude 
import obspy
import os
import math
from math import log10,sqrt
import numpy as np
from statistics import mean
from obspy import read,UTCDateTime
from obspy.geodetics import locations2degrees

ddirwaveform = '/Users/miao/Desktop/LOCFLOW/Data/waveform_sac'
stationdir = '../Data/station.dat'
phasedir = '../REAL/hypophase.dat' # hypoDD phase format
#phasedir = '../hypoDD/hypoDD.pha' # hypoDD phase format
catmag = './catalog_mag.txt'
g = open(catmag,'w')

# https://docs.obspy.org/_modules/obspy/signal/invsim.html
paz_wa = {'poles': [-6.283 + 4.7124j, -6.283 - 4.7124j],
                'zeros': [0 + 0j], 'gain': 1.0, 'sensitivity': 2080}

i=0
mags = []
with open(phasedir) as p:
    for lines in p:
        string = lines.split()[0]
        if string == '#':
            string, year, mon, day, hour, minute, tmp, lat, lon, dep, mag0, jk1, jk2, jk3, num = lines.split()
            sec,msec=tmp.split('.')
            date = year+mon+day
            if i>0:
                net_mag = np.median(mags)
                net_var = np.var(mags)
                g.write('{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(int(num)-1, year0, mon0, day0, hour0, minute0, tmp0, lat0, lon0, dep0, round(float(net_mag),2),round(float(net_var),2)))
                print('{} {} {} {} {} {} {} {} {} {} {} {}\n'.format(int(num)-1, year0, mon0, day0, hour0, minute0, tmp0, lat0, lon0, dep0, round(float(net_mag),2),round(float(net_var),2)))
                mags = []
            year0 = year
            mon0 = mon
            day0 = day
            hour0 = hour
            minute0 = minute
            tmp0 = tmp
            lat0 = lat
            lon0 = lon
            dep0 = dep
            i = 0
        else:
            station, tpick, tweight, phase = lines.split()
            with open(stationdir, "r") as f:
                for stations in f:
                    stlo, stla, net, sta, chan, ele = stations.split()
                    dist0 = 111.19*locations2degrees(float(lat), float(lon), float(stla), float(stlo))
                    if sta==station:
                        dist = math.sqrt(dist0**2+float(dep)**2)
                        if phase == 'P':
                            # cover the P phase and 3 sec after the S phase, just an example, change as needed.
                            tb = UTCDateTime(int(year),int( mon), int(day), int(hour), int(minute), int(sec), int(msec))+ float(tpick) - 0.5
                            te = tb + 0.73*(dist/6.0) + 3 #0.73*dist/6.0 is approximated to ts - tp
                        elif phase == 'S':
                            tb = UTCDateTime(int(year),int( mon), int(day), int(hour), int(minute), int(sec), int(msec)) + float(tpick) - 0.73*(dist/6.0) - 0.5
                            te = tb + 0.732*(dist/6.0) + 3
                        chann = chan[:2]+"[N,2]"
                        chane = chan[:2]+"[E,1]"
                        wavee = ddirwaveform+'/'+date+'/'+net+'.'+station+'.'+chane
                        waven = ddirwaveform+'/'+date+'/'+net+'.'+station+'.'+chann
                        tre = read(wavee,starttime=tb,endtime=te)
                        trn = read(waven,starttime=tb,endtime=te)
                        tre.simulate(paz_remove = None, paz_simulate = paz_wa)
                        trn.simulate(paz_remove = None, paz_simulate = paz_wa)
                        datatre = tre[0].data
                        datatrn = trn[0].data
                        #1000 is from meter to millimeter (mm) see Hutton and Boore (1987)
                        amp = (np.max(np.abs(datatre)) + np.max(np.abs(datatrn)))/2*1000
                        i = 1
                        #see Hutton and Boore (1987)
                        ml = log10(amp) + 1.110*log10(dist/100) + 0.00189*(dist-100)+3.0
                        #Just an example! please change into your magnitude formula.
                        mags.append(ml)
                        break
