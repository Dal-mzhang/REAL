#!usr/bin/perl -w
$mag = "2"; 
$num = 0;

$dir = "../REAL";

$together = "phase_allday.txt";
#Actually, it is the hypoDD input format	
open(OT,">$together");
	$file = "$dir/phase_sel.txt";
	open(JK,"<$file");  
	@par = <JK>;
	close(JK);
	foreach $file(@par){
		($test,$jk) = split(' ',$file);
		if($test =~ /^\d+$/){
			($jk,$year1,$mon1,$dd1,$time,$ot,$std,$lat,$lon,$dep) = split(' ',,$file);
			($hour,$min,$sec) = split('\:',$time);
			$num++;
			print OT "# $year1  $mon1  $dd1   $hour    $min    $sec    $lat    $lon    $dep     $mag     0.0     0.0    0.0   $num\n";
		}else{
			($net,$station,$phase,$traveltime,$pick,$amplitude) = split(' ',$file);
			print OT "$station $pick 1 $phase\n";
		}
	}
	
close(OT);
