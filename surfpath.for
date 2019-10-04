c--form inversion matrix-----------------------------------
	subroutine jcbm(nout,ipar,nor,qual,ndr,eqk,ndeq,sta,ndst,
     1    nrow,nz,tm3) 
	include "surfpath.par"
	include "datatype.inc"
	dimension val(300),nor(2,nry),tm3(3,nry),qual(nry)
	integer col(300)
	dimension nobk(nbk),azmbk(nbk),segbk(nbk),rcov(4,nbk),jazbk(4,nbk)
	type (station) sta(nst)
	type (earthquake) eqk(neq)
	type (partom) ipar

	ip=ipar%ipar;iwgh=ipar%iwgh; isplit=ipar%isplit; nr1gp=ipar%nr1gp
	nzvl=0.0;  nnzr=0;  nzcl=0;  bkw=0.0; iwdat=ipar%iwdat 
	npst=0
	do i=1,ndst
		npst=npst+sta(i)%nbin	
	end do
	nr=0;	nz=0
	rmh=6371.0-par0(1)
	pp=3.1415926/180.0
	vc=par0(2)

	nr=0;
	wdat=1.0;
	
	do ir=1,ndr
		if (iwdat.eq.0) then
			wdat=qual(ir);
		endif
		nr=nr+1;
		ie=nor(1,nr);
		is=nor(2,nr);
		elo=eqk(ie)%elo;ela=eqk(ie)%ela;edp=eqk(ie)%edp
		slo=sta(is)%slo;sla=sta(is)%sla;shg=sta(is)%shg
		dst=sphdstgc(elo,ela,slo,sla)*re
		az=sphazmgc(elo,ela,slo,sla)
		call bkofray(area,nx,ny,elo,ela,slo,sla,
     +			nbr,nobk,segbk,azmbk,re)
c		write(*,*)"ir=",ir," nbr=",nbr," nr=",nr," loc=",elo,ela,slo,sla
		if(nbr.le.0)then
			nr=nr-1
			goto 100
		end if
		tcal=ttime(ipar%ipar,modva,nx,ny,nbr,nobk,segbk,azmbk)
		tm3(2,nr)=tcal
		call raycov(nbk,nbr,nobk,segbk,azmbk,wdat,rcov,jazbk)
		call steq_pidx(ie,is,sta,ndst,eqk,ndeq,ips,ipe,nr1gp,isplit)
		call rwgh(dst,nbk,npar,npst,ipe,ips,nbr,nobk,segbk,bkw,ip)
		call jcbr(nbr,nbk,nobk,segbk,azmbk,npst,ipe,ips,rmh,
     +							nzr,col,val,ip)
		nnzr(nr)=nzr
		do iz=1,nzr
			nzcl(nz+iz)=col(iz)
			nzvl(nz+iz)=val(iz)*wdat;
		end do
		nz=nz+nzr				
100	end do
	nrow=nr
	do ibk=1,nbk
		bkazidx(ibk)=bk_azidx(jazbk(1,ibk))
	end do

	call datout2(nout,nx,ny,rcov)
	end
c------Get the index of the st/eq for one ray in the parameter list-------
	subroutine steq_pidx(ie,is,sta,ndst,eqk,ndeq,ipst,ipeq,n1gp,kdvd)
	include "datatype.inc"
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)

	if (kdvd.eq.1)then  !devided the # of rays by ray-number
		np0=0
		do i=1,is-1
			np0=np0+sta(i)%nbin						
		end do
		sta(is)%ntmp1=sta(is)%ntmp1+1
		if(sta(is)%ntmp1.gt.n1gp)then
			sta(is)%ntmp2=sta(is)%ntmp2+1			
			sta(is)%ntmp1=1			
		end if
		ipst=sta(is)%ntmp2+1
		if(ipst.gt.sta(is)%nbin)ipst=sta(is)%nbin
		ipst=ipst+np0

		np0=0
		do i=1,ie-1
			np0=np0+eqk(i)%nbin	
		end do
		eqk(ie)%ntmp1=eqk(ie)%ntmp1+1
		if(eqk(ie)%ntmp1.gt.n1gp)then
			eqk(ie)%ntmp2=eqk(ie)%ntmp2+1			
			eqk(ie)%ntmp1=1			
		end if
		ipeq=eqk(ie)%ntmp2+1
		if(ipeq.gt.eqk(ie)%nbin)ipeq=eqk(ie)%nbin
		ipeq=ipeq+np0
		return
	end if

	pp=3.1415926/180.0
	slo=sta(is)%slo;sla=sta(is)%sla;elo=eqk(ie)%elo;ela=eqk(ie)%ela
	sazm=sphazm_east(slo,sla,elo,ela)
	eazm=sphazm_east(elo,ela,slo,sla)

	np0=0
	do i=1,is-1
		np0=np0+sta(i)%nbin						
	end do
	do i=1,sta(is)%nbin
		if(sazm.ge.sta(is)%bin(1,i).and.sazm.le.sta(is)%bin(2,i))then
			ipst=np0+i
			goto 110
		end if			
	end do

110	np0=0
	do i=1,ie-1
		np0=np0+eqk(i)%nbin	
	end do
	do i=1,eqk(ie)%nbin
		if(eazm.ge.eqk(ie)%bin(1,i).and.eazm.le.eqk(ie)%bin(2,i))then
			ipeq=np0+i
			return
		end if			
	end do
	end 
c-----------form the weighting factor for cell and station delays-----------------
	subroutine formwgh(nbk,npar,sta,ndst,eqk,ndeq,bkw,wgh,ipar)
	include "datatype.inc"
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)
	type (partom)::ipar

	dimension bkw(npar),wgh(4)
	!-------calculate the weighting factor for Velocity+anisotropy--------
	if(ipar%iwgh.ne.0)then
		do ip=1,npar
			bkw(ip)=1.0
		end do
		return
	end if
	np=1
	if(ipar%ipar.eq.0)np=3
	do ibk=1,nbk
		ic=(ibk-1)*np		
		ww=sqrt(bkw(ic+1))
		if(ww.ne.0)then
			bkw(ic+1)=1.0/ww*wgh(1)
		end if
		if(ipar%ipar.eq.0)then
			ww=sqrt(bkw(ic+2))
			if(ww.ne.0)then
				bkw(ic+2)=1.0/ww*wgh(2)
			end if
			ww=sqrt(bkw(ic+3))
			if(ww.ne.0)then
				bkw(ic+3)=1.0/ww*wgh(2)
			end if
		end if
	end do
	!-------station and eq weighting not necessary for surface wave inversion
	return
	!-------calculate the weighting factor for station+eq--------
	ic=np*nbk
	do ist=1,ndst
		do j=1,sta(ist)%nbin
			ic=ic+1
			bb=bkw(ic)
			if(ipar%isplit.eq.0)then
				nrb=sta(ist)%nrb(j)
			elseif(ipar%isplit.eq.1)then
				nrb=ipar%nr1gp
				if(j.eq.sta(ist)%nbin.and.j.ne.1)then
					nrb=nrb+mod(sta(ist)%nrd,ipar%nr1gp)
				elseif(j.eq.sta(ist)%nbin.and.j.eq.1)then
					nrb=sta(ist)%nrd
				end if
			end if
			if(bb.ne.nrb)then
				write(*,*)"sorry,bb sould be equal to sta(ist).nrb(j)"
				!stop				
			end if
			if(bb.ne.0.0)bb=1.0/bb
			bkw(ic)=sqrt(bb)*wgh(3)
		end do
	end do
	do ieq=1,ndeq
		!ic=np*nbk+ndst+ieq
		do j=1,eqk(ieq)%nbin		
			ic=ic+1
			bb=bkw(ic)
			if(bb.ne.0.0)bb=1.0/bb
			bkw(ic)=sqrt(bb)*wgh(4)
		end do
	end do
	
	end
c--------calculating the weights for velocity,anisotropy----------------------------
	subroutine rwgh(dst,nbk,npar,npst,ipe,ips,nbr,nobk,segbk,bkw,ipar)
	dimension bkw(npar),nobk(nbk),segbk(nbk)
	np=1
	if(ipar.eq.0)np=3
	do ibr=1,nbr
		ib=nobk(ibr)
		ic=np*(ib-1)
		bkw(ic+1)=bkw(ic+1)+segbk(ibr)*dst !the weights for the velocity
		if(ipar.eq.0)then                   !the weights for the anisotropy
			bkw(ic+2)=bkw(ic+1)
			bkw(ic+3)=bkw(ic+1)
		end if
	end do
	!the weights for stations
c	bkw(nbk*np+ips)=bkw(nbk*np+ips)+1.0       
	!the weights for earthquakes
c	bkw(nbk*np+npst+ipe)=bkw(nbk*np+npst+ipe)+1.0  
	end
c------------------------------------
	subroutine raycovold(az,nbk,nbr,nobk,segbk,wdat,rcov)
	dimension nobk(nbk),segbk(nbk),rcov(4,nbk)
	do ibr=1,nbr
		ib=nobk(ibr)
		rcov(1,ib)=rcov(1,ib)+segbk(ibr)*wdat
		rcov(2,ib)=rcov(2,ib)+segbk(ibr)*sin(az)*wdat
		rcov(3,ib)=rcov(3,ib)+segbk(ibr)*cos(az)*wdat
		rcov(4,ib)=rcov(4,ib)+1.0*wdat
	end do
	end
c------------------------------------
	subroutine raycov(nbk,nbr,nobk,segbk,azmbk,wdat,rcov,jazm)
	dimension nobk(nbk),segbk(nbk),azmbk(nbk),rcov(4,nbk),jazm(4,nbk)
	pp=3.1415926/180.0;
	do ibr=1,nbr
		ib=nobk(ibr)
		az=azmbk(ibr)
		rcov(1,ib)=rcov(1,ib)+segbk(ibr)*wdat
		rcov(2,ib)=rcov(2,ib)+segbk(ibr)*sin(az)*wdat
		rcov(3,ib)=rcov(3,ib)+segbk(ibr)*cos(az)*wdat
		rcov(4,ib)=rcov(4,ib)+1.0*wdat

		iaz=int(az/pp/45.0)+1
		if (iaz.gt.8)then
			iaz=4
		elseif(iaz.gt.4)then
			iaz=iaz-4
		elseif(iaz.lt.1)then
			iaz=1
		endif
		jazm(iaz,ib)=jazm(iaz,ib)+1
	end do
	end
c----------------------------------------ccccccccccccccccccccccc
	subroutine jcbr(nbr,nbk,nobk,segbk,azmbk,ndst,ieq,ist,rmh,
     +							nzr,nzcol,val,ipar)
	dimension wgh(4),nzcol(nbr*3+2),val(nbr*3+2)
	dimension nobk(nbr),segbk(nbr),azmbk(nbr)

	inz=0
	if(ipar.eq.0)then  !velo+aniso+sta+eq
		nn=3
		do ib=1,nbr
			icol=(nobk(ib)-1)*nn
			bklen=segbk(ib)
			if(bklen.eq.0)goto 888
			inz=inz+1; 
			nzcol(inz)=icol+1;	val(inz)=bklen !velocity
			inz=inz+1; 
			nzcol(inz)=icol+2;	val(inz)=bklen*cos(2*azmbk(ib)) !Anisotropy A
			inz=inz+1; 
			nzcol(inz)=icol+3;	val(inz)=bklen*sin(2*azmbk(ib)) !Anisotropy B
888		end do
		! station parameter
c		inz=inz+1
c		icol=nbk*nn+ist
c		nzcol(inz)=icol;	val(inz)=1.0
		! source parameter
c		inz=inz+1
c		icol=nbk*nn+ndst+ieq
c		nzcol(inz)=icol;	val(inz)=1.0
	elseif(ipar.eq.1)then  !velo      +sta+eq
		nn=1
		do ib=1,nbr
			icol=(nobk(ib)-1)*nn
			bklen=segbk(ib)
			if(bklen.eq.0)goto 889
			inz=inz+1; 
			nzcol(inz)=icol+1;	val(inz)=bklen !velocity
889		end do
		! station parameter
c		inz=inz+1
c		icol=nbk*nn+ist
c		nzcol(inz)=icol;	val(inz)=1.0
		! source parameter
c		inz=inz+1
c		icol=nbk*nn+ndst+ieq
c		nzcol(inz)=icol;	val(inz)=1.0
	elseif(ipar.eq.2)then  
	!     aniso+sta+eq
	end if
	nzr=inz	
	end 
c-------------find the block in which one pn-ray pass-----------
	subroutine bkofray(area,nx,ny,eo,ea,so,sa,
     +		nbr,nobk,segbk,azmbk,re)
	dimension area(3,2),nobk(nx*ny),azmbk(nx*ny),segbk(nx*ny)
	real raypt(2,nx*ny), lo1,lo2,la1,la2
	nbk=nx*ny
	if(eo.lt.so)then
		lo1=eo; la1=ea
		lo2=so; la2=sa
	else
		lo1=so; la1=sa
		lo2=eo; la2=ea
	end if

c	write(*,*)lo1,la1,lo2,la2
	call raypts(area,nbk,lo1,la1,lo2,la2,npr,raypt)
	nbr=npr-1
	do ip=1,nbr
		p1o=raypt(1,ip);p1a=raypt(2,ip);
		p2o=raypt(1,ip+1);p2a=raypt(2,ip+1)
		!if(ip.le.3)write(*,*)p1o,p1a,p2o,p2a
		ibk=Nofbk(area,p1o,p1a,p2o,p2a);	
		if(ibk.eq.0)then
			nbr=0
			return
		end if
		nobk(ip)=ibk
		d2=seglength(area,raypt(1,ip),nx,ny,re)
		segbk(ip)=d2
		azmbk(ip)=sphazmgc(p1o,p1a,p2o,p2a)
	end do
	end
c!----the segment in every cell along the dipping moho------------------------
	function seglength(area,pt2,nx,ny,re)
	dimension pt2(2,2),area(3,2),dpt(2)
	dx=area(3,1)
	dy=area(3,2)
	pp=180.0/3.1415926
	ds=sphdstgc(pt2(1,1),pt2(2,1),pt2(1,2),pt2(2,2))
	r1=re; r2=re
	dhz=ds*re
	seglength=dhz;
	end
!-----depth interpolation with 3 points on one bigcicle and the depths of two-ebdpoints given------
	function p2depth(p1o,p1a,p1d,p2o,p2a,p2d,pmo,pma,re,nn)
	pi=3.1415926
	ds1=sphdst(pmo,pma,p1o,p1a)
	ds2=sphdst(pmo,pma,p2o,p2a)
	ds=sphdst(p1o,p2a,p2o,p2a)
	ds=ds1+ds2

	if(nn.eq.1)then
		w1=1.0/(ds1+0.000001)
		w2=1.0/(ds2+0.000001)
		p2depth=(p1d*w1+p2d*w2)/(w1+w2)
		return
	end if
	ac=ds

	sa=re-p1d
	sb=re-p2d
	sc=sa*sa+sb*sb-2.0*sa*sb*cos(ac)
	sc=sqrt(sc)

	sinaa=sa*sin(ac)/sc
	if(sinaa.gt.1.0)sinaa=1.0
	if(sinaa.lt.-1.0)sinaa=-1.0
	b2=3.1415926-ds2-asin(sinaa)

	rm=sinaa*sb/sin(b2)
	p2depth=re-rm
	end
!-----get the block number of one block-------------------
	function Nofbk(area,pt1o,pt1a,pt2o,pt2a)
	dimension area(3,2)
	nx=(area(2,1)-area(1,1))/area(3,1)
	ny=(area(2,2)-area(1,2))/area(3,2)
	pta=(pt1a+pt2a)/2.0;pto=(pt1o+pt2o)/2.0
	ix=int((pto-area(1,1))/area(3,1))
	iy=int((pta-area(1,2))/area(3,2))
	Nofbk=iy*nx+ix
	if(nofbk.gt.nx*ny)nofbk=0
	end 
!----find cross points of great circle of longitude and latitude lines-------
	subroutine raypts(area,nbk,lo1,la1,lo2,la2,np,raypt)
	dimension area(3,2)
	real lo1,lo2,la1,la2, raypt(2,nbk), path(2,nbk),trgl(3,2)
	!---when azimuth is equal to zero000000000000000
	if(abs(lo1-lo2).le.0.01)then
		ny1=int((la1-area(1,2))/area(3,2))
		ny2=int((la2-area(1,2))/area(3,2))
	   	if(ny1.lt.ny2)then
			np=1;	raypt(1,np)=lo1;	raypt(2,np)=la1;
			do ia=ny1+1,ny2
				np=np+1
				ptla=ia*area(3,2)+area(1,2)
				raypt(1,np)=lo1;	raypt(2,np)=ptla;
			end do
			np=np+1
			raypt(1,np)=lo1;	raypt(2,np)=la2;
		else
			np=1;	raypt(1,np)=lo2;	raypt(2,np)=la2;
			do ia=ny2+1,ny1
				np=np+1
				ptla=ia*area(3,2)+area(1,2)
				raypt(1,np)=lo2;	raypt(2,np)=ptla;
			end do
			np=np+1
			raypt(1,np)=lo2;	raypt(2,np)=la1;
		end if
		return
	end if
	!---when azimuth is not zero----------------
	nx0=int((lo1-area(1,1))/area(3,1))
	nx1=int((lo2-area(1,1))/area(3,1))

	pi=3.1415926
	pp=3.1415926/180.0
	xlo1=lo1*pp
	yla1=pi/2-la1*pp
	xlo2=lo2*pp
	yla2=pi/2-la2*pp
	azm=sphazm(lo1,la1,lo2,la2)
	np=1
	path=0.0
	path(1,np)=lo1; path(2,np)=la1; 
	do ix=nx0+1,nx1  
	!find the cross points of the great circle 
	!and all longitudes between two end-points
		ptlo=area(1,1)+ix*area(3,1)
		xptlo=ptlo*pp
c		ptla=gcptla(lo1,la1,lo2,la2,ptlo)
		ang1=azm; side=yla1; ang2=xptlo-xlo1
		gcp=sphside(ang1,side,ang2)
		ptla=(pi/2-gcp)/pp

		if(inarea(area,ptlo,ptla).eq.0)then
			np=0
			return
		end if
	    doa=abs(ptlo-path(1,np))+abs(ptla-path(2,np))
	    if(doa.gt.0.000001)then
			np=np+1
			path(1,np)=ptlo
			path(2,np)=ptla
		end if
	end do
	if(path(1,np).ne.lo2)then
		np=np+1
		path(1,np)=lo2;	path(2,np)=la2
	end if
	np0=np
	np=1
	!find the cross points of the great circle 
	!and all latitudes between two end-points and rearrange
	raypt=0.0
	raypt(1,1)=path(1,1);raypt(2,1)=path(2,1)
	do ip=1,np0-1
		p1o=path(1,ip)
		p1a=path(2,ip)
		p2o=path(1,ip+1)
		p2a=path(2,ip+1)
		ny0=int((p1a-area(1,2))/area(3,2))
		ny1=int((p2a-area(1,2))/area(3,2))
		if(ny1.gt.ny0)then
			do iy=ny0+1,ny1
				ptla=area(1,2)+iy*area(3,2)
				if(ptla.ne.p2a)then
					yptla=pi/2-ptla*pp
					ptlo=gcptlo(p1o,p1a,p2o,p2a,ptla)
					
					if(inarea(area,ptlo,ptla).eq.0)then
						np=0
						return
					end if
				    doa=abs(ptlo-raypt(1,np))+abs(ptla-raypt(2,np))
				    if(doa.gt.0.000001)then
						np=np+1
						raypt(1,np)=ptlo
						raypt(2,np)=ptla
					end if
				end if
			end do
		elseif(ny1.lt.ny0)then
			do iy=ny0,ny1+1,-1
				ptla=area(1,2)+iy*area(3,2)
				if(ptla.ne.p1a)then
					yptla=pi/2-ptla*pp
					ptlo=gcptlo(p1o,p1a,p2o,p2a,ptla)

					if(inarea(area,ptlo,ptla).eq.0)then
						np=0
						return								 
					end if
				    doa=abs(ptlo-raypt(1,np))+abs(ptla-raypt(2,np))
					if(doa.gt.0.000001)then
						np=np+1
						raypt(1,np)=ptlo
						raypt(2,np)=ptla
					end if
				end if
			end do
		end if
	    doa=abs(p2o-raypt(1,np))+abs(p2a-raypt(2,np))
		if(doa.gt.0.000001)then
			np=np+1
			raypt(1,np)=p2o
			raypt(2,np)=p2a
		end if
	end	 do
	end
!------add the regularization to jacobe matrix------------------
	subroutine rgljcb(ibk,np,ip,wei,nz)
	include "surfpath.par"
	nz=nz+1
	icol=(ibk-1)*np+ip
	nzcl(nz)=icol
	nzvl(nz)=-1.0*wei
	end
