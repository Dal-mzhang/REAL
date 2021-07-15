#!/bin/bash
phasein="../REAL/phase_sel.txt"
stationin="../Data/station.dat"
velocityin="../REAL/tt_db/itvel.nd"

####step 1##### create the phase file
python mk_inputfile.py $phasein $stationin > hypoinput.arc

####step 2##### create the velocity model, pay more attentation
#python mk_velmodel.py $velocityin
#please run mk_velmodel.py to convert the Vp and Vs models
#then mannually adjust them
#Note: the hypoinverse cannot allow two layers have the same velocity
#If you get constant depths, you have this issue! 
#Slightly ajdust the model (e.g., increase vel a little bit with depth)

####step 3####### run hypoinverse
hyp1.40 <hyp.command

###convert to readable format
cat hypoOut.arc|gawk '{if(substr($0,1,4) == 2016)print substr($1,1,4),substr($1,5,2),substr($1,7,2),substr($1,9,2),substr($1,11,2),substr($1,13,4)/100,substr($1,17,2)+substr($0,20,4)/6000,substr($0,24,3)+substr($0,28,4)/6000,substr($0,33,4)/100}'|gawk '{printf "%4s%2s%2s %2d %2d %5.2f %7.4f %8.4f %5.2f %d\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,NR}'>new.cat
