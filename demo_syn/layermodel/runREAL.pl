#!/usr/bin/perl -w
$day = "14";

$D = "2016/10/$day";
$R = "0.2/20/0.02/2/5";
$G = "1.4/20/0.01/1";
$V = "6.2/3.3";
$S = "5/0/18/1/0.5/0.0/1.0";

$dir = "./pk";
$station = "../station.dat";
$ttime = "../tt_db/ttdb.txt";

system("time REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime");
print"REAL -D$D -R$R -G$G -S$S -V$V $station $dir $ttime\n";

#A few of picks are lost (i.e., < 120) because we only keep the most reliable pick within a time window.
#For real case, REAL will select picks based on their weighting factors
