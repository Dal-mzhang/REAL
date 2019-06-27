#!/usr/bin/perl -w
$phasein = "../VELEST/phase_allday.txt";
$event = "../VELEST/new.cat";

$phaseout = "italy.pha"; #phase format for hypoDD

open(JK,"<$phasein");
@par = <JK>;
close(JK);

open(EV,"<$event");
@eves = <EV>;
close(EV);

open(EV,">$phaseout");
foreach $file(@par){
    chomp($file);
    ($test,$jk) = split(' ',$file);
    if($test eq "#"){
		($jk,$year,$month,$day,$hour,$min,$sec,$lat,$lon,$dep,$jk,$jk,$jk,$jk,$num) = split(' ',,$file);
			$out = 0;
		foreach $eve(@eves){
			chomp($eve);
			($date,$hh,$mm,$ss,$lat1,$lon1,$dep1,$num0,$gap,$res) = split(" ",$eve);
			$num1 = $num0*100;
            #combine original picks &  improved initial locations that obtained with VELEST
			if(abs($num-$num1)<1 && abs($hour*3600 + $min*60 + $sec - $hh*3600 - $mm*60 - $ss) < 1.5){
				print EV "# $year $month $day $hour $min $sec $lat1 $lon1 $dep1 0 0 0 0 $num\n";
				$out = 1;
			}
		}
	}else{
		if($out>0){
			printf EV "$file\n";
		}
    }
}
close(EV);
