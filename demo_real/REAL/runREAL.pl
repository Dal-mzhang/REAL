#!/usr/bin/perl -w
#input example for REAL 1.3
#should be changed for your case 
$year = "2016";
$mon = "10";
$day = "14";

# -D(nyear/nmon/nday/lat_center)
# latitude center is required to make lat. and lon. consistent in km
$D = "$year/$mon/$day/42.75";
# -R(rx/rh/tdx/tdh/tint[/gap/GCarc0/latref0/lonref0]])
# station gap constraint is strongly recommended when station coverage is poor
$R = "0.2/20/0.02/2/3/250/2";
# -G(trx/trh/tdx/tdh)
$G = "1.4/20/0.01/1";
# -V(vp0/vs0/[s_vp0/s_vs0/ielev])
$V = "6.2/3.3";
# -S(np0/ns0/nps0/npsboth0/std0/dtps/nrt[/rsel/ires])
$S = "3/2/5/2/0.5/0.2/1.5";

$dir = "../Pick/$year$mon$day";
$station = "../Data/station.dat";
$ttime = "./tt_db/ttdb.txt";

system("REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";
