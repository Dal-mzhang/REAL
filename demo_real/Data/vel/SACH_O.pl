#!/usr/bin/env perl
use warnings;
#Here we define 00h00m00.00s of one day as our reference time and mark the origin time as ZERO.

for($day=14;$day<=14;$day++){
$year = "2016";
$mon = "10";
$dir1 = "$year$mon$day";
$dir2 = "$year$mon$day";

$jday = 288+($day-14);
$hour = "00";
$min = "00";
$sec0 = "00";
$msec = "00";

@sac = glob "$dir1//*";
open(SAC,"|sac> jk");
for($i=0;$i<@sac;$i++){
	chomp($sac[$i]);
	($jk,$name) = split("//",$sac[$i]);
	print SAC "r $sac[$i]\n";
	print SAC "ch o gmt $year $jday $hour $min $sec0 $msec\n";
	print SAC "evaluate to tmp 0 - &1,o\n";
	print SAC "ch allt %tmp% iztype io\n";
	print SAC "ch o 0\n";
	print SAC "w $dir2/$name\n";
}
print SAC "q\n";
close(SAC);

open(SAC,"|sac > jk");
print SAC "cut 0 86400\n";
print SAC "r $dir2/*\n";
print SAC "w over\nq\n";
close(SAC);
unlink "jk";

}
