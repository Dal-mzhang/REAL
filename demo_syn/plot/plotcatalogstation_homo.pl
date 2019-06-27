#!/usr/bin/perl -w
$event = "../layermodel/catalog_use.txt";
$locate2 = "../homomodel/catalog_sel.txt";
$sta = "../station.dat";
$R = "12.5/13.95/42.22/43.28";
$J = "M6.5i";
$ps = "Italy_homo.ps";
$cpt = "globe";
$evlo = 13.25;
$evla = 42.75;

$location = "location";
$station = "station";
$relocate2 = "locate2";
#process 
open(JK,"<$event");
@par = <JK>;
close(JK);

open(FK,"> temp");
foreach $_(@par){
        chomp($_);
		($ot,$lon,$lat,$depth,$mag) = split(" ",$_);
        ($date,$time) = split("T",$ot);
        ($year,$month,$day) = split('\-',$date);
        ($hour,$min,$sec) = split('\:',$time);
		print FK "$lon $lat $depth $mag\n";
}
close(FK);
`more temp | sort -k3r > $location`;

open(JK,"<$sta");
@par = <JK>;
close(JK);
open(FK,"> $station");
foreach $_(@par){
        chomp($_);
		($lon,$lat,$net,$sta0,$comp,$elev) = split(" ",$_);
		print FK "$lon $lat $net.$sta0\n";
}
close(FK);

open(JK,"<$locate2");
@par = <JK>;
close(JK);
open(FK,"> $relocate2");
foreach $_(@par){
        chomp($_);
		($num,$yy,$mm,$dd,$otime,$time,$std,$lat,$lon,$dep,$mag,$jk) = split(" ",$_);
		print FK "$lon $lat\n";
}
close(FK);

#globe setting
`gmt gmtset FORMAT_GEO_MAP DDD:mm:ssF`;
`gmt gmtset FONT_ANNOT_PRIMARY 18p`;

#start plot
`gmt psbasemap -P -R$R -J$J -Ba0.3  -BWesN -K > $ps`;
`gmt pscoast -R -J -Dh -W1/0.8p -I1/0.5p -Na/0.5p -K -O >> $ps`;


`more $location | gmt psxy -R -J -Sa0.28c -Gblue -K -O >> $ps `;
`more $relocate2 | gmt psxy -R -J -Sa0.28c -Gred -K -O >> $ps `;

#plot stations
`gmt psxy $station -R -J -St0.5c -Gblack -K -O >> $ps `;

#plot circle
$temp = $evla-0.47; #in degree
`echo $evlo $evla 100 | gmt psxy -J -R -SE- -K -O -W2p,gray >> $ps`;
`echo $evlo $temp 50 km | gmt pstext -J -R -F+f20p,black -K -O >> $ps`;

#plot whole map
$R = "6.5/19/36.5/45.5";
$J = "M1.6i";
$B = "a1000/a1000";
`gmt gmtset MAP_FRAME_TYPE plain`;
`gmt psbasemap -P -R$R -J$J -B$B -K -O -X0.1i -Y0.1i  >> $ps`;
`gmt pscoast -R -J -Dh -W1/0.5p -I1/0.5p -K -O >> $ps`;
open(GMT,"|gmt psxy -R -J -W2p,black -K -O >> $ps");
print GMT "12.5 42.22\n";
print GMT "12.5 43.28\n";
print GMT "13.95 43.28\n";
print GMT "13.95 42.22\n";
print GMT "12.5 42.22\n";
close(GMT);


`gmt psxy -R -J -T -O >> $ps`;
`rm gmt.*  temp $location $relocate2 $station`;
