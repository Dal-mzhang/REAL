#!/usr/bin/perl -w
$day = "14";

$D = "2016/10/$day";
$R = "0.1/20/0.02/2/5";
$G = "1.4/20/0.01/1";
$V = "6.2/3.3";
$S = "3/2/5/2/0.5/0.0/1.5";

$dir = "./pk";
$station = "../station.dat";
$ttime = "../tt_db/ttdb.txt";

system("time ../../bin/REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"../../bin/REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";
