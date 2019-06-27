#!/usr/bin/perl -w
$year = "2016";
$mon = "10";
$day = "14";

$D = "$year/$mon/$day";
$R = "0.2/20/0.02/2/5";
$G = "1.4/20/0.01/1";
$V = "6.2/3.3";
$S = "5/0/18/1/0.5/0.5/1.3/1.8";

$dir = "../Pick/$year$mon$day";
$station = "../Data/station.dat";
$ttime = "./tt_db/ttdb.txt";

system("REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";
