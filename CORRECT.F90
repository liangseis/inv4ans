	!-------correction-------------------------------------
	subroutine correct(eq,st,slw,v0,dth,dte)
	dimension eq(3),st(3),ph(23),dh(23)
	data ph /0.0,10.11,16.12,20.42,24.29,27.53,31.00,33.58,36.48,39.33,42.15,44.55,47.37, &
		50.17,53.01,55.50,58.46,61.51,65.11,68.52,73.11,78.50,90.00/
	data dh /7,6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-14/
	!-------station altitude correction-------
	aa=180*slw*v0/(3.1415926*6371.0)
	pp=3.1415926/180.0
	aa=aa**2
	aa=1.0-aa
	if(aa.ge.0.0)then
		aa=sqrt(aa)
	else
		aa=0.0
	end if
	aa=aa*st(3)/1000.0
	dth=aa/v0
	!----------------------------------------
	re=eq(1)/180.0*3.1415926
	pe=3.1415926/2.0-eq(2)*pp
	rs=st(1)/180.0*3.1415926
	ps=3.1415926/2.0-st(2)*pp

	elo=eq(1)
	ela=eq(2)
	slo=st(1)
	sla=st(2)
	dlt=sphdst(elo,ela,slo,sla)
	dlt=dlt/pp

	if(dlt.le.9)then
		fd=0.01
	elseif(dlt.le.13)then
		fd=0.02
	elseif(dlt.le.22)then
		fd=0.03
	elseif(dlt.le.41)then
		fd=0.04
	elseif(dlt.le.63)then
		fd=0.05
	elseif(dlt.le.72)then
		fd=0.06
	elseif(dlt.le.83)then
		fd=0.07
	else
		fd=0.08
	end if

	nph=23
	yla=eq(2); ieq=isearch(yla,nph,ph); heq=dh(ieq)
	yla=st(2); ist=isearch(yla,nph,ph); hst=dh(ist)
	dte=fd*(heq+hst)
	end
	!-----------------------------------------
	function isearch(yla,nph,ph)
	dimension ph(nph)
	yla=abs(yla)
	do ip=1,nph-1
		if(yla.ge.ph(ip).and.yla.le.ph(ip+1))then
			isearch=ip
			return
		end if
	end do
	isearch=nph
	end
