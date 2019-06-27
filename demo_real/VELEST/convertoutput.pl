#!/usr/bin/perl -w
#161015  0 6 49.51 42.8676N  13.0841E   4.14   0.00    102      0.24
$relocate = "final.CNV";
$relo = "new.cat";
$dele = "delet.cat";

open(JK,"<$relocate");
@par = <JK>;
close(JK);

open(OT,">$relo");
open(DE,">$dele");
foreach $_(@par){
	chomp($_);
	if(substr($_,0,4) eq "1610"){
	$date = substr($_,0,6);
	$hour = substr($_,7,2);
	$min = substr($_,9,2);
	$sec = substr($_,12,5);
	$lat = substr($_,18,7);
	$lon = substr($_,27,8);
	$dep = substr($_,37,6);
	$mag = substr($_,45,5);
	$az = substr($_,54,3);
	$res = substr($_,62,5);
    if(substr($_,25,1) eq 'S'){$lat = -1*$lat;}
    if(substr($_,35,1) eq 'W'){$lon = -1*$lon;}
	if($az <= 200 && $res < 0.6){
		print OT "$date $hour $min $sec $lat $lon $dep $mag $az $res\n";
	}else{
		print DE "$date $hour $min $sec $lat $lon $dep $mag $az $res\n";
	}
	}
}
close(OT);
close(DE);
