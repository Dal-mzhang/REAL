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

system("REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";
