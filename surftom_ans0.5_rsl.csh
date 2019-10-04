#! /bin/csh -f

set fl="surftom.ini"
set ptho="./"
set pthi="../../dspt2/"

set ptn=3.0
#foreach prd ('050' '060' '070' '080' '090' '100' '110' '120' '130' '140' '150' '160' '170' '180' '190' '200' '210' '220' '230' '240' '250' '260' '270' '280' '290' '300' '310' '320' '330' '340' '350' '360' '370' '380' '390' '400' '410' '420' '430' '440' '450' '460' '470' '480' '490' '500')
foreach prd ('100' '200' '300' '400' '450' '500')
#foreach prd ('370' '380')
	set flo = surftom.$prd.ini
	echo "Group Velocity inversion of etibet:" >! $flo
	echo '"'$ptho'"'>> $flo
	echo '"'${pthi}tgrp_p${prd}etibet.txt'"' >> $flo
	if ($prd <= '1900') then
		echo "1.0 $ptn 0.3" >> $flo      ## //0--real inversion only; other--real+resol test//details see notes #1
	endif	
        echo "0" >> $flo 		## //0--one time inversion //other--bootstrap (no. of iterations).
	echo "0 1000" >> $flo 		## //1=split the rays of one station into several groups.0=no//
	echo "1 1 5 1 1" >> $flo 	## //number of bins for stations/earthquak delays
	echo "0" >> $flo		## //icrr---station and earthquake corrections  0==no correction  b4 inversion.
	echo "0" >> $flo		## //parameter included 0-with anisotropy/1-without anisotropy/
	echo "0" >> $flo		## //ds/v/dv inversion  0--dt=x*ds   1--absoluate velocity(not used?)  2--dt=x*dv/(-v0*c0);
	echo "0  0" >> $flo		## //iwgh--weighting cell or not? iwdat--weighting data or not 1--no weighting
	echo "0 0.0 0.0" >> $flo	## //1=stability test,two runs; 0-no test but can add error: cutoff, deviation 
	if ($prd <= '200') then
               set sm=0.5
		set sm2=1.0
        else if ($prd <= '250') then
		set sm=1.0
		set sm2=2.0
        else if ($prd <= '300') then
                set sm=1.5
		set sm2=3.0
        else if ($prd <= '350') then
                set sm=2.0
		set sm2=4.0
        else if ($prd < '400') then
                set sm=3.0
		set sm2=6.0
        else if ($prd < '450') then
                set sm=4.0
		set sm2=8.0
        else
                set sm=5.0
		set sm2=9.9
        endif
	echo $sm $sm2 >> $flo

	echo "15.0 15.0 15.0 15.0" >> $flo	## //weights for Velocity/anisotropy/stations/earthquake
	echo "1 -75 75"	>> $flo		## //kdt dtmin dtmax : in seconds: kdt==1, limit applied, ==0: limit not applied
	echo "1 1 35"	>> $flo		## //kds dsmin, dsmax: in degrees: kds==1, limit applied, ==0: limit not applied
	echo "200" >> $flo		##  //maximum iteration times
	echo "rlomod8_ptn${ptn}box_p${prd}.dat" >> $flo
	echo "rlooth8_ptn${ptn}box_p${prd}.dat" >> $flo
	echo "/The inverted model file: Used to check real inversion results only " >> $flo
	echo "rlomod_p150s10.dat" >> $flo
	echo "/Spike Check Points " >> $flo
	echo "0 -0.3 0.0" >> $flo
	cp $flo $fl
	surftom2_etp0.50
end
