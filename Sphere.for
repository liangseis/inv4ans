c-----find the latitude by a given longitude on a big-circle
	function gcptla(lo1,la1,lo2,la2,ptlo)
	!all inputs are in degrees
	real la1,lo1,la2,lo2,dist

	if(lo2.lt.lo1)then
		xx=lo2; lo2=lo1; lo1=xx
		yy=la2; la2=la1; la1=yy
	end if

	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=lo1*pp
	yla1=pi/2-la1*pp
	xlo2=lo2*pp
	yla2=pi/2-la2*pp
	xptlo=ptlo*pp

	dist=sphdst(lo1,la1,lo2,la2) !inputs are in degrees
	side1=yla1; side2=dist; side3=yla2
	ang=sphang12(side1,side2,side3)

	ang1=ang; side=yla1; ang2=xptlo-xlo1
	gcptla=sphside(ang1,side,ang2)
	gcptla=(pi/2-gcptla)/pp
	end
c-----find the longitude by a given latitude on a big-circle
	function gcptlo(lo1,la1,lo2,la2,ptla)
	real la1,lo1,la2,lo2,dist

	if(lo2.lt.lo1)then
		xx=lo2; lo2=lo1; lo1=xx
		yy=la2; la2=la1; la1=yy
	end if
	if(abs(ptla-la1).le.0.0001)then
		gcptlo=lo1; 
		return
	elseif(abs(ptla-la2).le.0.0001)then
		gcptlo=lo2; 
		return
	end if
	pi=3.1415926
	pp=pi/180.0
	xlo1=lo1*pp
	yla1=pi/2-la1*pp
	xlo2=lo2*pp
	yla2=pi/2-la2*pp
	yptla=pi/2-ptla*pp

	dist=sphdst(lo1,la1,lo2,la2) !inputs are in degrees
	side1=yla1; side2=dist; side3=yla2
	ang=sphang12(side1,side2,side3)

	side1=yla1; side2=yptla; ang2=ang
	dlo=xlo2-xlo1
	ang3=sphang3(side1,side2,ang2,dlo)
	gcptlo=(ang3+xlo1)/pp
	end
c--------------calculate the epicenter distnce--------
	function sphdstgc(elo,ela,slo,sla)
	!gc---Convert GeoGraphic to GeoCentric
	!find the spherical distance between two points
	!all inputs are in degree
	call gg2gc(ela,elac)
	call gg2gc(sla,slac)
	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=elo*pp
	yla1=pi/2-elac*pp
	xlo2=slo*pp
	yla2=pi/2-slac*pp

	ang=abs(xlo2-xlo1)
	if(ang.gt.pi)ang=ang-pi
	
      dst=cos(yla1)*cos(yla2)+sin(yla1)*sin(yla2)*cos(ang)
	if(dst.lt.-1.0)dst=-1.0
	if(dst.gt.1.0)dst=1.0
	dst=acos(dst)
	sphdstgc=dst
	end
c--------------calculate the azimuth angle of two points on sphere--------
	function sphazmgc(xo1,ya1,xo2,ya2)
	!gc---Convert GeoGraphic to GeoCentric
	!all inputs are in degree, output are in radian degree

	do0=abs(xo1-xo2)*3.1415926/180.0;
	if(do0.le.0.001)then
		sphazmgc=0.0
		return
	end if

	call gg2gc(ya1,ya1c)
	call gg2gc(ya2,ya2c)
	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=xo1*pp
	yla1=pi/2-ya1c*pp
	xlo2=xo2*pp
	yla2=pi/2-ya2c*pp
	ang=abs(xlo2-xlo1)
	if(ang.le.0.00001)then
	   cosx=1.0
	 else
	    cosx=cos(ang)
	 end if
	if(ang.gt.pi)ang=ang-pi
      dst=cos(yla1)*cos(yla2)+sin(yla1)*sin(yla2)*cosx
	if(dst.le.0.00001)dst=0.00001;
	dst=acos(dst)

	tmp=sin(yla1)*sin(dst)
	cosaz=(cos(yla2)-cos(yla1)*cos(dst))/tmp
      sinaz=sin(yla2)*sin(xlo2-xlo1)/sin(dst)
      azimuth=atan2(sinaz,cosaz)
      if(azimuth.lt.0.0) azimuth =azimuth+2*pi
      sphazmgc=azimuth

c	write(*,*)"SPHAZM",yla1,yla2,ang,dst,azimuth
	end
c------------------------------------------
	subroutine gg2gc(grla,cela)
		b2asq=0.9933055
		radperdeg=3.1415926/180.0
		cela=atan(b2asq*tan(radperdeg*grla))/radperdeg
	end
c--------------calculate the epicenter distnce--------
	function sphdst(elo,ela,slo,sla)
	!find the spherical distance between two points
	!all inputs are in degree
	pi=3.1415926
	pp=3.1415926/180
	xlo1=elo*pp
	yla1=pi/2-ela*pp
	xlo2=slo*pp
	yla2=pi/2-sla*pp

	ang=abs(xlo2-xlo1)
	if(ang.gt.pi)ang=ang-pi
	
      dst=cos(yla1)*cos(yla2)+sin(yla1)*sin(yla2)*cos(ang)
	dst=acos(dst)

	sphdst=dst
	end
c--------------calculate the azimuth angle of two points on sphere--------
	function sphazm_east(xo1,ya1,xo2,ya2)
	!find the spherical distance between two points
	!all inputs are in degree, output are in radian degree
	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=xo1*pp
	yla1=pi/2-ya1*pp
	xlo2=xo2*pp
	yla2=pi/2-ya2*pp

	if(abs(xo1-xo2).le.0.001)then
		sphazm=0.0
		return
	end if
	ang=abs(xlo2-xlo1)
	if(ang.gt.pi)ang=ang-pi
	
      dst=cos(yla1)*cos(yla2)+sin(yla1)*sin(yla2)*cos(ang)
	dst=acos(dst)

	if(dst.eq.0.0)write(*,*)xo1,ya1,xo2,ya2,dst,sin(dst)
	tmp=sin(yla1)*sin(dst)
	cosaz=(cos(yla2)-cos(yla1)*cos(dst))/tmp
      sinaz=sin(yla2)*sin(xlo2-xlo1)/sin(dst)
      azimuth=atan2(sinaz,cosaz)
      if(azimuth.lt.0.0) azimuth =azimuth+2*pi
      sphazm_east=azimuth/pp-90.0
	if(sphazm_east.lt.0.0)sphazm_east=sphazm_east+360.0
	end
c--------------calculate the azimuth angle of two points on sphere--------
	function sphazm(xo1,ya1,xo2,ya2)
	!find the spherical distance between two points
	!all inputs are in degree, output are in radian degree
	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=xo1*pp
	yla1=pi/2-ya1*pp
	xlo2=xo2*pp
	yla2=pi/2-ya2*pp

	if(abs(xo1-xo2).le.0.001)then
		sphazm=0.0
		return
	end if
	ang=abs(xlo2-xlo1)
	if(ang.gt.pi)ang=ang-pi
	
      dst=cos(yla1)*cos(yla2)+sin(yla1)*sin(yla2)*cos(ang)
	dst=acos(dst)

	if(dst.eq.0.0)write(*,*)xo1,ya1,xo2,ya2,dst,sin(dst)
	tmp=sin(yla1)*sin(dst)
	cosaz=(cos(yla2)-cos(yla1)*cos(dst))/tmp
      sinaz=sin(yla2)*sin(xlo2-xlo1)/sin(dst)
      azimuth=atan2(sinaz,cosaz)
      if(azimuth.lt.0.0) azimuth =azimuth+2*pi
      sphazm=azimuth
	end

c---find the angle between the first two sides with 3 sides given----
	function sphang12(side1,side2,side3)
	p=(side1+side2+side3)/2
	cosa2=sqrt(sin(p)*sin(p-side3)/sin(side1)/sin(side2))
	sphang12=2*acos(cosa2)
	end
c--find angle with 2 sides and one angle given-----------------
	function sphang3(side1,side2,ang2,dlo)
	!find the spherical angle between side1 and side2
	!all inputs are in radian degree
	iflag=0
	pi=3.1415926
	sina1=sin(side1)*sin(ang2)/sin(side2)
	ang1=asin(sina1)
22	ctanAB=cos((ang1+ang2)/2)/sin((ang1+ang2)/2)
	tan=cos((side2-side1)/2)*ctanAB/cos((side2+side1)/2)
	ang3=atan(tan)
	if(ang3.lt.0)ang3=ang3+pi
	ang3=ang3*2
	if (ang3.gt.dlo)then
		if (iflag.eq.1)then
			sphang3=ang3
			return
		else
			iflag=1
			ang1=pi-ang1
			goto 22
		end if
	end if
	sphang3=ang3
	end
c--find one side with 2 angles and their in-between side given---------------
	function sphside(ang1,side,ang2)
	!find the spherical subtense of angle1
	!all inputs are in radian degree
	ss0=cos((ang2-ang1)/2)/(cos((ang2+ang1)/2))*tan(side/2)
	ss0=atan(ss0)
	ss1=sin((ang2-ang1)/2)/(sin((ang2+ang1)/2))*tan(side/2)
	ss1=atan(ss1)
	sphside=ss0-ss1
	if(sphside.lt.0)sphside=sphside+3.1415926
	end
c--find the angle of the first vertex of a triangle with 3 vertex points given---------------
	function sphang13p(px1,py1,px2,py2,px3,py3)
	!find the spherical subtense of angle1 all inputs are in radian degree
	!all inputs are in degree
	!output in degree
	side12=sphdst(px1,py1,px2,py2)
	side13=sphdst(px1,py1,px3,py3)
	side23=sphdst(px2,py2,px3,py3)
	sphang13p=sphang12(side1,side2,side3)*180.0/3.1415926
	end
!-----one point is in one spherical rectangular area??----------------------------
	function inrect(rec,pnt)
	dimension rec(2,4),pnt(2)
	real ya(2)
	inrect=0
	iya=0; ya(1)=0.0; ya(2)=0.0
	if(pnt(1).ge.rec(1,1).and.pnt(1).le.rec(1,4))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,1),rec(2,1),rec(1,4),rec(2,4),pnt(1))
	elseif(pnt(1).ge.rec(1,4).and.pnt(1).le.rec(1,1))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,4),rec(2,4),rec(1,1),rec(2,1),pnt(1))
	end if

	if(pnt(1).ge.rec(1,2).and.pnt(1).le.rec(1,3))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,2),rec(2,2),rec(1,3),rec(2,3),pnt(1))
	elseif(pnt(1).ge.rec(1,3).and.pnt(1).le.rec(1,2))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,3),rec(2,3),rec(1,2),rec(2,2),pnt(1))
	end if
	if(iya.ge.2)goto 111
	if(pnt(1).ge.rec(1,1).and.pnt(1).le.rec(1,2))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,1),rec(2,1),rec(1,2),rec(2,2),pnt(1))
	elseif(pnt(1).ge.rec(1,2).and.pnt(1).le.rec(1,1))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,2),rec(2,2),rec(1,1),rec(2,1),pnt(1))
	end if
	if(iya.ge.2)goto 111
	if(pnt(1).ge.rec(1,3).and.pnt(1).le.rec(1,4))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,3),rec(2,3),rec(1,4),rec(2,4),pnt(1))
	elseif(pnt(1).ge.rec(1,4).and.pnt(1).le.rec(1,3))then
		iya=iya+1
		ya(iya)=gcptla(rec(1,4),rec(2,4),rec(1,3),rec(2,3),pnt(1))
	end if

	if(iya.lt.2)return
111	if(ya(1).le.ya(2))then
		if(pnt(2).ge.ya(1).and.pnt(2).le.ya(2))inrect=1
	else
		if(pnt(2).ge.ya(2).and.pnt(2).le.ya(1))inrect=1
	end if
	end
c-----------------------------------------------------------
	function sphaexc(px1,py1,px2,py2,px3,py3)
c	!for the given 3 sperical points, find the angle-excess
c	!delta=sum of the three spherical angles minus the pi
c	!in/output are in degree
	sd12=sphdst(px1,py1,px2,py2)
	sd13=sphdst(px1,py1,px3,py3)
	sd23=sphdst(px2,py2,px3,py3)
	ang1=sphang12(sd23,sd13,sd12)*180.0/3.1415926
	ang2=sphang12(sd13,sd12,sd23)*180.0/3.1415926
	ang3=sphang12(sd12,sd13,sd23)*180.0/3.1415926
	delta=ang1+ang2+ang3-180.0
	sphaexc=0.0;
	end

	subroutine sph2flt(vs,zs,vf,zf,re)
		rr=(re-zs)/re
		zf=-1.0*re*log(rr)
		vf=vs/rr
	end

	subroutine flt2sph(vf,zf,vs,zs,re)
		rs=re*exp(-1.0*zf/re)
		zs=re-rs
		vf=vs*rs/re
	end
