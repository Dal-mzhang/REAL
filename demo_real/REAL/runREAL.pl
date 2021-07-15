#!/usr/bin/perl -w
#should be changed for your case 
$year = "2016";
$mon = "10";
$day = "14";

# -D(nyear/nmon/nday/lat_center)
# latitude center is required to make lat. and lon. consistent in km
$D = "$year/$mon/$day/42.75";
# -R(rx/rh/tdx/tdh/tint[/gap/GCarc0/latref0/lonref0]])
#grid size determine the speed, use different grid size
#-S should be different
#$R = "0.1/20/0.02/2/5"; # example: use a small grid size
$R = "0.1/20/0.04/2/5"; # example: use a large grid size
# -G(trx/trh/tdx/tdh)
$G = "1.4/20/0.01/1";
# -V(vp0/vs0/[s_vp0/s_vs0/ielev])
$V = "6.2/3.3";
# -S(np0/ns0/nps0/npsboth0/std0/dtps/nrt[drt/nxd/rsel/ires])
#$S = "3/2/10/3/0.5/0.2/1.6/0.5/0.2/4"; # example: use a small grid size 
$S = "3/2/12/3/0.5/0.2/1/0.25/0.2/4"; # example: use a large grid size

$dir = "../Pick/$year$mon$day";
$station = "../Data/station.dat";
$ttime = "./tt_db/ttdb.txt";

system("../../bin/REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"../../bin/REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";

#how to evaluate your parameters:
#1. check P and S t_dist curve (if too narrow, increase nrt; vice versa)
#2. make sure you lose very few events after the second selection (if lose too many, increase thresholds, e.g., np0, ns0, nps0 and/or npsboth0)
#3. go to event_verify and check those worst events (large station gap and traveltime residual) and make sure they are ture events, 
#   otherwise, increase your thresholds.
#4. go to the ../hypoinverse dir, run hypoinverse, unstable events usually show REMARK "-", use their ID to check their waveform in event_verify
#   if they are false events, increase your thresholds.
    
#tips:
#1. decrease drt (0 is the best) will improve the stability but will cost more time, 0 means no initiating picks will be removed.
#   It is time afforable (a couple of minutes) if your total number of grids is only about 1000 (i.e.,10x10x10).
#2. nrt and grid size (tdx and tdh) trade off, when you use small grid size, choose a large nrt (but may < 2.0)
#   similarly, when you use large grid size, please use a small nrt. You don't want to remove too many associated picks in the initiating pick pool
#   and increase the risk. Thus, please decrease the drt as well.
#3. nxd and rsel can be used to remove those suspicious events and picks, respectively.
#4. you may get less events when you use loose thresholds compared to strict ones when the drt > 0.
#5. the main purpose is assocaition, not try to refine the grid to have a better location.
#6. more picks or stations, more strict thresholds.
#7. common problems: too slow -> the total number of girds is too big
#                    flase events -> too loose thresholds
#                    bad performance -> didn't use optimal parameters for different grid size
