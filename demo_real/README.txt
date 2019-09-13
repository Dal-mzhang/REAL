Steps 1&2 are time-consuming. I would suggest you start with step 3 by assuming you have needed picks.

1. download and process seismic data (skip it to save time)
cd Data
cd vel
python ingv_waveform_IV_vel.py #download seismic data from INGV
python ingv_waveform_YR_vel.py #download seismic data from IRIS
perl SAC_O.pl #change 00:00:00 as origin time
cd ../
perl transferWA.pl #create waveforms with corrected amplitude to calculate local magnitude (Wood-Anderson)
(It is optional. You don't need to do this if you don't want to estimate magnitude)

2. pick the first arrivals of P and S phases using STA/LTA method (skip it to save time)
cd Pick (unzip 20161014.zip and move to step 3 to save time)
python trigger_p_amp.py #create P picks and amplitudes using STA/LTA
python trigger_s_amp.py #create S picks and amplitudes using STA/LTA
(I have been working on earthquake association and location using machine-learning picks and REAL)
(Compared with STA/LTA, both detection ability and locaiton accruacy are significantly improved)
(But machine-learning picker is not available in this package)

3. Associate and locate events using STA/LTA picks and REAL
cd tt_db 
python taup_tt.py # bulid a traveltime table
(You may change the searching range and interval. Remember to change them in runREAL.pl as well)
(You may change the velocity model following TauP model format)
perl runREAL.pl
(Amplitudes in pick files are used for magnitude estimation but they are optional. If you don't have amplitudes, set as 0.0)
cd t_dist
awk -f pha_t-dist.awk #create travel time vs. distance
t_dist.m (matlab) #roughly estimate association result from travel time vs. distance

4. refine earthquake locations using VELEST
cd VELEST
mergetogether_phase.pl #collect all phase inforamtion and initial locations determined by REAL
perl converformat.pl #convert format (see velest manual)
velest (type the command, you may change velest.cmn file in your case)
perl convertoutput.pl #convert location format. Further select events based on station gap and residual 
(some false detections may occur around or close to the boundary of the study area, remove events with large station gaps!!)
plot_3dscatter.m (matlab) #plot location distribution in 3-D

5. further improve earthquake locations using hypoDD
cd hypoDD
bash createstation.sh #create station file for hypoDD
perl selectphase.pl #create phase file for hypoDD (combine original picks and improved locations by VELEST)
ph2dt ph2dt.inp #create paired traveltime difference
hypoDD hypoDD.inp #you may change parameter settings for your case. You may use the improved velocity model that obtained with VELEST (velout.mod).
plot_3dscatter.m (matlab) #plot locations in 3-D

