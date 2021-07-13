#!/usr/bin/perl -w
$station = "../Data/station.dat";

$stause = "italy.sta"; #station format for VELEST
open(JK,"<$station");
@par = <JK>;
close(JK);

$p1=0;
$p2=1;
$p3=1;
$v1=0.00;
$v2=0.00;
open(OT,">$stause");
#print OT "(a4,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2)\n"; #old code
print OT "(a6,f7.4,a1,1x,f8.4,a1,1x,i4,1x,i1,1x,i3,1x,f5.2,2x,f5.2)\n"; #new version
foreach $_(@par){
	chomp($_);
	($lon,$lat,$net,$sta,$comp,$elev) = split(" ",$_);
	#if(length($sta)>4){$sta = substr($sta,1,4);} 
	$vsn = "N";$vew = "E";
	if($lat < 0.0){$vsn = "S";$lat = -1*$lat;}
	if($lon < 0.0){$vew = "W";$lon = -1*$lon;}
	printf OT "%-6s%7.4f%1s %8.4f%s %4d %1d %3d %5.2f  %5.2f\n",$sta,$lat,$vsn,$lon,$vew,$p1,$p2,$p3,$v1,$v2;
	$p3++;
}
print OT "\n";
close(OT);


##################################################
#Note: VELEST doesn't have event index
#Here I use magnitude to represent its index
#The index would be used in hypoDD procedure as well
#Hope I didn't confuse you
#################################################
$phasein = "./phase_allday.txt";
$phaseout = "italy.pha"; # phase format for VELEST ISED=1
$phasecat = "italy.cat"; # catalog format for VELEST

open(JK,"<$phasein");
@par = <JK>;
close(JK);

open(EV,">$phaseout");
open(CT,">$phasecat");
foreach $file(@par){
    chomp($file);
    ($test,$jk) = split(' ',$file);
    if($test eq "#"){
		($jk,$year,$month,$day,$hour,$min,$sec,$lat,$lon,$dep,$jk,$jk,$jk,$jk,$num) = split(' ',,$file);
		$year = substr($year,2,2); # VELEST format 
		$vsn = "N";$vew = "E";
		if($lat < 0.0){$vsn = "S"; $lat = -1*$lat;} # VELEST format
		if($lon < 0.0){$vew = "W"; $lon = -1*$lon;}
		
		$mag = $num/100;
        	#(3i2,1x,2i2,1x,f5.2,1x,f7.4,a1,1x,f8.4,a1,1x,f7.2,2x,f5.2)
		print EV "\n";
		printf EV "%2d%2d%2d %2d%2d %5.2f %7.4f%s %8.4f%s %7.2f  %5.2f\n",$year,$month,$day,$hour,$min,$sec,$lat,$vsn,$lon,$vew,$dep,$mag;
		printf CT "%2d%2d%2d %2d%2d %5.2f %7.4f%s %8.4f%s %7.2f  %5.2f\n",$year,$month,$day,$hour,$min,$sec,$lat,$vsn,$lon,$vew,$dep,$mag;
	}else{
        ($station,$tpick,$jk,$phase) = split(' ',$file);
		$iwt = "0";
		#if(length($station)>4){$station = substr($station,1,4);} #again, old format...
        #(2x,a6,2x,a1,3x,i1,3x,f6.2)
        printf EV "  %-6s  %-1s   %1d   %6.2f\n",$station,$phase,$iwt,$tpick;
    }
}
print EV "\n";
close(EV);
close(CT);
