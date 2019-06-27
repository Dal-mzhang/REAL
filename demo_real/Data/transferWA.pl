#!/usr/bin/perl -w

$year = "2016";
$mon = "10";
$dayb = "14";
$daye = "14";

$wa = "./wa";
if(!-e $wa){
	`mkdir $wa`;
}
for($day=$dayb;$day<=$daye;$day++){
	$dir = "$year$mon$day";
	@sacs = glob "./vel/$dir/*.SAC";
		
	if(!-e "$wa/$dir"){
		mkdir "$wa/$dir";	
	}

	open(SAC,"|sac > jk");
	foreach $_(@sacs){
		chomp($_);
		print SAC "r $_\n";
		#print SAC "rmean\nrtr\ntaper width 0.0001\n"; #taper has been applied in previous data processing
		print SAC "bp c 0.05 15 n 4 p 2\n";
		print SAC "transfer from vel to general n 2 f 0.8 d 0.7 m 2080\n"; #transfer to Wood-Anderson instrument
		print SAC "w $_.wa\n";
	}
	print SAC "q\n";
	close(SAC);
	`mv ./vel/$dir/*.wa $wa/$dir/`;
    unlink "jk";
}
