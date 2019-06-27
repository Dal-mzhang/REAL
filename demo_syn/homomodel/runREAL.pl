#!/usr/bin/perl -w
$day = "14";

$D = "2016/10/$day";
$R = "0.2/20/0.02/2/5";
$V = "6.2/3.3";
$S = "5/0/18/1/0.5/0.0/1.0";

$dir = "../layermodel/pk";
$station = "../station.dat";
$ttime = "../tt_db/ttdb.txt";

system("time REAL -D$D -R$R -S$S -V$V $station $dir");
print"REAL -D$D -R$R -S$S -V$V $station $dir\n";

#Note that no traveltime table and -G 
#A few of picks are lost (i.e., < 120) because we only keep the most reliable pick within a time window.
#For real case, REAL will select picks based on their weighting factors
