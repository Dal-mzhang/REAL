Steps 1&2 are time-consuming. I would suggest you start with step 3 by assuming you have needed picks.

1. download and process seismic data (skip it to save time)
cd Data
cd vel
python ingv_waveform_IV_vel.py #download seismic data from INGV
python iris_waveform_YR_vel.py #download seismic data from IRIS
perl SAC_O.pl #change 00:00:00 as origin time
cd ../
perl transferWA.pl #create waveforms with corrected amplitude to calculate local magnitude (Wood-Anderson)
(It is optional. You don't need to do this if you don't want to estimate magnitude)

2. pick the first arrivals of P and S phases using STA/LTA method (skip it to save time)
cd Pick (unzip 20161014.zip and move to step 3 to save time)
python trigger_p_amp.py #create P picks and amplitudes using STA/LTA
python trigger_s_amp.py #create S picks and amplitudes using STA/LTA
(see Pick_PhaseNet, it shows how to get machine-learning picks)
(Compared with STA/LTA, both detection ability and locaiton accruacy are significantly improved)
(PhaseNet can be found at https://github.com/wayneweiqiang/PhaseNet)

3. Associate and locate events using STA/LTA picks and REAL
cd REAL
cd tt_db 
python taup_tt.py # bulid a traveltime table
(You may change the search range and interval. Remember to change them in runREAL.pl as well)
(You may change the velocity model following TauP model format)
cd ../
perl runREAL.pl
(Amplitudes in pick files are used for magnitude estimation but they are optional. If you don't have amplitudes, set as 0.0)
cd t_dist
awk -f pha_t-dist.awk #create travel time vs. distance
t_dist.m (matlab) #roughly estimate association result from travel time vs. distance
cd event_verify # verify events with the smallest number of picks using python scripts to evaluate your parameters
new version REAL generates simulated annealing locations, which can be directly used in hypoDD.

4. refine earthquake locations using VELEST
cd VELEST
perl mergetogether_phase.pl #collect all phase inforamtion and initial locations determined by REAL
perl converformat.pl #convert format (see velest manual)
(please prepare your own velocity model x.mod following the instruction in velest manual)
velest (type the command to run velest, you may need to change parameters in velest.cmn file in your case such as reference point lat and lon)
(note that the reference lontitude "olon" is reversed, e.g., 13.25E should be -13.25 here, see velest manual)
(there two options in velest: 1) "isingle" = 1, you will update location alone. recommend ittmax=99 (a large number) and invertratio=0; It is fine if you see the code stops with errors in the end; 2) "isingle" = 0, you will update location, velocity, and station correction. recommend ittmax=9 and invertratio=3. remember to have consistent "neqs" with the number of events in this case, here see demo velest.cmn_0, move velest.cmn_0 to velest.cmn before you run the code)
perl convertoutput.pl #convert location format and further select events based on station gap and residual 
(some false detections may occur around or close to the boundary of the study area, remove events with large station gaps or large residuals)
plot_3dscatter.m (matlab) #plot location distribution in 3-D


5. refine earthquake locations using hypoinverse
cd hypoinverse
bash run_hypoinverse.sh #format conversion scripts are provided

6. further improve earthquake locations using hypoDD
you may use the REAL's SA locations, VELEST locations, or hypoinverse locations as initial locations.
here use velest as example
cd hypoDD
bash createstation.sh #create station file for hypoDD
perl selectphase.pl #create phase file for hypoDD (combine original picks and improved locations by VELEST)
ph2dt ph2dt.inp #create paired traveltime difference
hypoDD hypoDD.inp #you may change parameters in your case. You may use the improved velocity model (isingle=0) that obtained with VELEST (velout.mod).
plot_3dscatter.m (matlab) #plot locations in 3-D

7. calculate magnitude
cd Magnitude
python calc_mag.py #an example for magnitude estimation (assume you have waveforms after removing response)
