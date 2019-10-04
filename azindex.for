	!-----------calculating the azimuth distribution level of every block----------
	function bk_azidx(jaz)
	dimension jaz(4)
	nnn=0;
	ann=0.0;
	do i=1,4
		ann=ann+jaz(i)
		if (jaz(i).gt.0) then
			nnn=nnn+1
		end if
	end do
	if(ann.eq.0.0)then
		bk_azidx=3.0
		return
	end if

	az1=0;
	do i=1,2
		ff=(jaz(i)+jaz(i+2))
		if (ff.le.0) then
			az1=az1+1.0;
		else
			az1=az1+abs(jaz(i)-jaz(i+2))/ff
		end if
	end do

	bk_azidx=az1/nnn;
	end
	!-----------calculating the azimuth distribution level of stations----------
	function steq_azidx(jaz)
	dimension jaz(8)
	nnn=0;
	ann=0.0;
	do i=1,8
		ann=ann+jaz(i)
		if (jaz(i).gt.0) then
			nnn=nnn+1
		end if
	end do
	if(ann.eq.0.0.or.nnn.eq.0)then
		azidx=99.0
		return
	end if

		az1=0;  		!rays are evenly distributed in 8 az-bins, then az1=0
		do i=1,4
			ff=(jaz(i)+jaz(i+4))
			if (ff.le.0) then
				az1=az1+1.0;
			else
				az1=az1+abs(jaz(i)-jaz(i+4))/ff
			end if
		end do
		az1=az1/4

			az2=0;  !rays are evenly distributed in 8 az-bins, then az2=0
			do i=1,8
				j1=i+3; j2=i+5
				if(j2.gt.8) then
					j2=j2-8
				end if
				ff=(2*jaz(i)+jaz(j1)+jaz(j2))
				if (ff.le.0) then
					az2=az2+1.0;
				else
					az2=az2+abs(2*jaz(i)-jaz(j1)-jaz(j2))/ff
				end if
			end do
			az2=az2/8

			az3=0;  !rays are evenly distributed in 8 az-bins, then az3=0
			do i=1,8
				j1=i+2; j2=i+6
				if(j2.gt.8) then
					j2=j2-8
				end if
				ff=(2*jaz(i)+jaz(j1)+jaz(j2))
				if (ff.le.0) then
					az3=az3+1.0;
				else
					az3=az3+abs(2*jaz(i)-jaz(j1)-jaz(j2))/ff
				end if
			end do
			az3=az3/8

			az4=0;
			do i=1,8
				j1=i+1; j2=i+7
				if(j2.gt.8) then
					j2=j2-8
				end if
				ff=(2*jaz(i)+jaz(j1)+jaz(j2))
				if (ff.le.0) then
					az4=az4+1.0;
				else
					az4=az4+abs(2*jaz(i)-jaz(j1)-jaz(j2))/ff
				end if
			end do
			az4=az4/8


	if (nnn.gt.4)then
		azidx=(az1+az2+az3+az4)/4.0
	else
		azidx=(az1+az2+az3+az4)/4.0
	end if
	steq_azidx=azidx !/nnn
	end
