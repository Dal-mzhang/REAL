import matplotlib.pyplot as plt
import obspy
import os
from obspy import read,UTCDateTime
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, classic_sta_lta
import warnings
import numpy as np
from obspy import UTCDateTime

year = 2016
mon = 10
day = 14

date = str(year)+str(mon)+str(day)+'/'


def main():
    if not os.path.exists(date):
        os.makedirs(date)
    ddir = '../Data/vel/'
    ddirwa = '../Data/wa/'

    with open("../Data/station.dat", "r") as f:
        for station in f:
            stlo, stla, net, sta, chan, elev = station.split()
            chane = chan[:2]+"E"
            chann = chan[:2]+"N"

            wave_e = ddir+date+net+'.'+sta+'.'+chane+'.'+'SAC'
            wave_n = ddir+date+net+'.'+sta+'.'+chann+'.'+'SAC'
            wave_e_wa = ddirwa+date+net+'.'+sta+'.'+chane+'.'+'SAC.wa'
            wave_n_wa = ddirwa+date+net+'.'+sta+'.'+chann+'.'+'SAC.wa'
        
            try:
                ste = read(wave_e)
                stn = read(wave_n)
                stewa = read(wave_e_wa)
                stnwa = read(wave_n_wa)
            
                tre = ste[0]
                trn = stn[0]
            
                tre_wa = stewa[0]
                trn_wa = stnwa[0]
        
                tre.detrend('demean') 
                tre.detrend('linear') 
                trn.detrend('demean') 
                trn.detrend('linear') 
                tre.filter(type="bandpass",freqmin=2.0,freqmax=15.0,zerophase=True)
                trn.filter(type="bandpass",freqmin=2.0,freqmax=15.0,zerophase=True)
                df = tre.stats.sampling_rate
                tstart = tre.stats.starttime - UTCDateTime(year, mon, day, 0, 0, 0)
                #print(tstart);
                output = './'+date+net+'.'+sta+'.'+'S.txt'

                # Characteristic function and trigger onsets
                cft = recSTALTAPy_h(tre.data,trn.data, int(0.2 * df), int(2.5 * df))
                on_of = trigger_onset(cft, 4.0, 2.0)
            
                tre_wa.detrend('demean') 
                tre_wa.detrend('linear') 
                trn_wa.detrend('demean') 
                trn_wa.detrend('linear') 
                tre_wa.filter(type="bandpass",freqmin=0.2,freqmax=10.0,zerophase=True)
                trn_wa.filter(type="bandpass",freqmin=0.2,freqmax=10.0,zerophase=True)
                datatre = tre_wa.data
                datatrn = trn_wa.data

                # Output the triggered 
                f = open(output,'w')
                i = 0
                while(i<len(on_of)):
                    trig_on = on_of[i,0]
                    trig_of = on_of[i,1]
                    trig_off = int(trig_of + (trig_of - trig_on)*4.0)
                    amp = max(max(abs(datatre[trig_on:trig_off])),max(abs(datatrn[trig_on:trig_off])))*1000
                    if max(cft[trig_on:trig_of]) > 6.0:
                        f.write('{} {} {}\n'.format(tstart+trig_on/df,max(cft[trig_on:trig_of]),amp))
                    i=i+1
                f.close()
            except:
                print('no data in ',wave_e,wave_n)


def recSTALTAPy_h(a, b, nsta, nlta):
    """
    Recursive STA/LTA written in Python.

    .. note::

        There exists a faster version of this trigger wrapped in C
        called :func:`~obspy.signal.trigger.recSTALTA` in this module!

    :type a: NumPy ndarray
    :param a: Seismic Trace
    :type nsta: Int
    :param nsta: Length of short time average window in samples
    :type nlta: Int
    :param nlta: Length of long time average window in samples
    :rtype: NumPy ndarray
    :return: Characteristic function of recursive STA/LTA

    .. seealso:: [Withers1998]_ (p. 98) and [Trnkoczy2012]_
    """
    try:
        a = a.tolist()
    except:
        pass

    try:
        b = b.tolist()
    except:
        pass
    ndat = len(a)
    # compute the short time average (STA) and long time average (LTA)
    csta = 1. / nsta
    clta = 1. / nlta
    sta = 0.
    lta = 1e-99  # avoid zero devision
    charfct = [0.0] * len(a)
    icsta = 1 - csta
    iclta = 1 - clta
    for i in range(1, ndat):
        sq = a[i] ** 2 + b[i] ** 2
        sta = csta * sq + icsta * sta
        lta = clta * sq + iclta * lta
        charfct[i] = sta / lta
        if i < nlta:
            charfct[i] = 0.
    return np.array(charfct)

if __name__ == '__main__':
    main()
