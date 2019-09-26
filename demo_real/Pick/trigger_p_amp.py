import matplotlib.pyplot as plt
import obspy
import os
from obspy import read,UTCDateTime
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta

year = 2016
mon = 10
day = 14

date = str(year)+str(mon)+str(day)+'/'

if not os.path.exists(date):
    os.makedirs(date)

ddir = '../Data/vel/'
ddirwa = '../Data/wa/'

with open("../Data/station.dat", "r") as f:
    for station in f:
        stlo, stla, net, sta, chan, elev = station.split()
        chanz = chan[:2]+"Z"
        chann = chan[:2]+"N"
        chane = chan[:2]+"E"
        
        wave = ddir+date+net+'.'+sta+'.'+chanz+'.'+'SAC'
        wavee = ddirwa+date+net+'.'+sta+'.'+chane+'.'+'SAC.wa'
        waven = ddirwa+date+net+'.'+sta+'.'+chann+'.'+'SAC.wa'

        try:
            st = read(wave)
            ste = read(wavee)
            stn = read(waven)
            
            tr = st[0]  
            tr.detrend('demean') 
            tr.detrend('linear') 
            tr.filter(type="bandpass",freqmin=2.0,freqmax=24.0,zerophase=True)
            df = tr.stats.sampling_rate
            tstart = tr.stats.starttime - UTCDateTime(year, mon, day, 0, 0, 0)
            #print(tstart)
            output = './'+date+net+'.'+sta+'.'+'P.txt'

            # Characteristic function and trigger onsets, see ObsPy website
            cft = recursive_sta_lta(tr.data, int(0.1 * df), int(2.5 * df))
            on_of = trigger_onset(cft, 6.0, 2.0)

            # Corrected amplitude (local magnitude)
            tre = ste[0]
            tre.detrend('demean') 
            tre.detrend('linear') 
            tre.filter(type="bandpass",freqmin=0.2,freqmax=10.0,zerophase=True)
            datatre = tre.data
        

            trn = stn[0]
            trn.detrend('demean') 
            trn.detrend('linear') 
            trn.filter(type="bandpass",freqmin=0.2,freqmax=10.0,zerophase=True)
            datatrn = trn.data

            # Output the triggered 
            f = open(output,'w')
            i = 0
            while(i<len(on_of)):
                trig_on = on_of[i,0]
                trig_of = on_of[i,1]
                trig_off = int(trig_of + (trig_of - trig_on)*4.0)
                #1000 is from meter to millimeter (mm) see Hutton and Boore (1987)
                amp = max(max(abs(datatre[trig_on:trig_off])),max(abs(datatrn[trig_on:trig_off])))*1000
                if max(cft[trig_on:trig_of]) > 10.0:
                    f.write('{} {} {}\n'.format((tstart+trig_on/df),max(cft[trig_on:trig_of]),amp))
                i=i+1
            f.close()
        except:
            print('no data in',wave)
