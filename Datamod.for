c------read the records from rec.dat---------------------------
	subroutine surfdata(luinv,fmck,frc,nout,nor,tm3,qual,ndeq,ndst,
     +							ndr,eqk,sta,ipar)
	include "surfpath.par"
	include "datatype.inc"
	dimension nor(2,nry),tm3(3,nry),qual(nry)
	dimension nobk(nbk),segbk(nbk),azmbk(nbk)
	dimension ee(3),ss(3),cmode(6)
	character stn*5,eqn*5,frc*80,fmd0*80,str*6,fmck*80
	real mdck(nbk,3)
	dimension jazst(8,nst),jazeq(8,neq)
	type (station)::sta(nst)
	type (earthquake)::eqk(neq)
	type (partom)::ipar

	cmode=ipar%cmode; icrr=ipar%icrr
	nor=0; tm3=0.0; ir=0; ndeq=0; ndst=0; ndr=0;
	!--------1-----------read the model first
	iter=0
	open(1,file=frc,status="old")
99		read (1,"(a6)")str							  
		if(str.ne." /rays")goto 99
		read(1,*,end=886)par0(3),nall
		par0(1)=0.0
		!par0(1)=-re*log((re-par0(1))/re)
		par0(2)=6.3
		!par0(3)=par0(3)*(re/(re-par0(1)))

		CALL inimod(nx,ny,nbk,modva,par0)
		adp=par0(1)
		vc=par0(2)
		vmean=par0(3)
		write(*,*)par0(1),par0(2),par0(3)
		if(cmode(1).ne.0.and.cmode(1).ne.99.and.cmode(1).ne.7)then
			call ckmod(fmck,luinv,mdck,cmode)  !resolution
		elseif (cmode(1).eq.7)then
			call SpikeMod(fmck,luinv,mdck,cmode,ipar%subreg)
		end if
	!--------2-----------read the records---------------------------
		pp=3.1415926/180.0
		ia=0
		ir=0
100		ela=0.0
		read(1,*,end=886)
     +				ela,elo,edp,sla,slo,shg,dst,tim,quali,err,eqn,stn
c		read(1,*,end=886)
c	+				ela,elo,edp,sla,slo,shg,dst,tim,quali,eqn,stn
c		err=0.0;
		if(ela.eq.0.0)goto 886

c		if (ela.eq.29.95.and.elo.eq.-95.833) goto 100
c		----------------------------------
c		quali=quali*10.0/(10.0+err): take into account the errors in reading
c	    ----06/18/2006---C. Liang----------
		quali=quali*10.0/(10.0+err/2.0);
		ia=ia+1
c		if(mod(ia,10).eq.0)write(*,*)ia
		!--2.1--Choose data by research area/dist/residual/apparent velocity
		if(inarea(area,elo,ela).eq.0)goto 100
		if(inarea(area,slo,sla).eq.0)goto 100
		if(dst.lt.dsmin.or.dst.gt.dsmax)goto 100
		call bkofray(area,nx,ny,elo,ela,slo,sla,nbr,nobk,
     +							    segbk,azmbk,re)
		if(nbr.le.0)goto 100
		if(quali.lt.0.0)goto 100 
		tcal=ttime(ipar%ipar,modva,nx,ny,nbr,nobk,segbk,azmbk) 
c		call correct(ee,ss,slw,v0,dth,dte)
c		tim=tim-(dth+dte)
		dtt=tim-tcal

		write(*,*)" in surfdata: iray,tim,tcal,dtt=",ia,tim,tcal,dtt

		if(abs(dtt).ge.dtmax.and.kdt==1)goto 100
		if (abs(dtt)<err .and. err>7.5) goto 100

c		if(abs(dtt).ge.35 .and. dst<600.0)goto 100
	    if (abs(dtt).gt.15.0.and.abs(dtt).le.25.0)then	
c			quali=quali*0.85
		elseif (abs(dtt).gt.25.0.and.abs(dtt).le.35.0)then	
c			quali=quali*0.7
		elseif (abs(dtt).gt.35.0)then	
c			quali=quali*0.5
		endif
121		vpna=dst*111.199/tim

c		if(vpna.gt.vmean+0.5.or.vpna.lt.vmean-0.5)goto 100
		if(cmode(1).ne.0.and.cmode(1).ne.99)then  !checkboard only
			tck=ttime(ipar%ipar,mdck,nx,ny,nbr,nobk,segbk,azmbk) 
			tim=tck
			dtt=tim-tcal
		end if				

		slw=1.0/par0(3); v0=par0(2)
		ee(1)=elo;ee(2)=ela;ee(3)=edp;ss(1)=slo;ss(2)=sla;ss(3)=shg
		write(nout,*)tim,tcal,dtt
		!--2.2--Number of rays for every event and station
		dddt=dddt+abs(dtt)
		dt2=dt2+dtt*dtt
		ir=ir+1;
		tm3(1,ir)=tim; 	tm3(2,ir)=tcal;
c		if(icrr.eq.1)tm3(1,ir)=tim-ecrr-scrr;

		mbin0=ipar%mbin0
		da=360.0/mbin0
		jazms=int(sphazm_east(slo,sla,elo,ela)/da)+1
		if (jazms .lt.1) then
			jazms=1;
		elseif(jazms .gt.mbin0)then
			jazms=mbin0;
		end if
		jazme=int(sphazm_east(elo,ela,slo,sla)/da)+1
		if (jazme .lt.1) then
			jazme=1;
		elseif(jazme .gt.mbin0)then
			jazme=mbin0;
		end if
		if(ir.eq.1)then
			ndeq=ndeq+1; ndst=ndst+1
			eqk(1)%elo=elo;eqk(1)%ela=ela;eqk(1)%edp=edp
			sta(1)%slo=slo;sta(1)%sla=sla;sta(1)%shg=shg;
			eqk(1)%eqn=eqn; sta(1)%stn=stn;
			
c			sta(1).crr=scrr;eqk(1).crr=ecrr
			eqk(1)%nrd=eqk(1)%nrd+1; 	nor(1,ir)=1
			sta(1)%nrd=sta(1)%nrd+1;  nor(2,ir)=1
c			eqk(1).iyr=iyr; eqk(1).imn=imn; eqk(1).idy=idy; 
c			eqk(1).ihr=ihr; eqk(1).imi=imi; eqk(1).sec=sec; 
			eqk(1)%iehb=0; eqk(1)%iset=0; 
			qual(1)=quali; 
			jazst(jazms,ndst)=jazst(jazms,ndst)+1;
			jazeq(jazme,ndeq)=jazeq(jazme,ndeq)+1;
		else
			if(elo.ne.eqk(ndeq)%elo.or.ela.ne.eqk(ndeq)%ela)then
				ndeq=ndeq+1	
				eqk(ndeq)%crr=0.0
				eqk(ndeq)%elo=elo;eqk(ndeq)%ela=ela;eqk(ndeq)%edp=edp
c				eqk(ndeq).iyr=iyr;eqk(ndeq).imn=imn;eqk(ndeq).idy=idy
c				eqk(ndeq).ihr=ihr;eqk(ndeq).imi=imi;eqk(ndeq).sec=sec 
				eqk(ndeq)%iehb=0; eqk(ndeq)%iset=0;
				eqk(ndeq)%eqn=eqn;
			end if
			do is=1,ndst
				if(slo.eq.sta(is)%slo.and.sla.eq.sta(is)%sla)then
					nor(2,ir)=is
					sta(is)%nrd=sta(is)%nrd+1
					jazst(jazms,is)=jazst(jazms,is)+1;
					goto 885
				end if
			end do			
			ndst=ndst+1 
			nor(2,ir)=ndst
			sta(ndst)%nrd=sta(ndst)%nrd+1 	
			sta(ndst)%slo=slo	
			sta(ndst)%sla=sla 
			sta(ndst)%shg=shg
c			sta(ndst).crr=scrr; 
			sta(ndst)%stn=stn
			jazst(jazms,ndst)=jazst(jazms,ndst)+1
885			eqk(ndeq)%nrd=eqk(ndeq)%nrd+1
			jazeq(jazme,ndeq)=jazeq(jazme,ndeq)+1
			nor(1,ir)=ndeq
			qual(ir)=quali
		end if

		goto 100
886	close(1)
	!--3--Choose data by number of rays of each station and event
	ndr=ir
	
	!split rays into different bins by azimuth
	do is=1,ndst
		sta(is)%azidx=steq_azidx(jazst(1,is))
		mbin0=ipar%mbin0;mbst=ipar%mbst;minr=ipar%minrst;
		if(ipar%isplit.eq.0)then
		call bingrp8azm
     + 		(jazst(1,is),mbin0,mbst,minr,nbin,sta(is)%bin,sta(is)%nrb)
		sta(is)%nbin=nbin; 
		elseif(ipar%isplit.eq.1)then
			sta(is)%nbin=int(sta(is)%nrd/ipar%nr1gp);
			if(sta(is)%nbin.eq.0)sta(is)%nbin=1
		end if
	end do
	do ie=1,ndeq
		eqk(ie)%azidx=steq_azidx(jazeq(1,ie))
		mbin0=ipar%mbin0;mbeq=ipar%mbeq;minr=ipar%minreq;
		if(ipar%isplit.eq.0)then
		call bingrp8azm
     +		(jazeq(1,ie),mbin0,mbeq,minr,nbin,eqk(ie)%bin,eqk(ie)%nrb)
		eqk(ie)%nbin=nbin; 
		elseif(ipar%isplit.eq.1)then
			eqk(ie)%nbin=int(eqk(ie)%nrd/ipar%nr1gp);
			if(eqk(ie)%nbin.eq.0)eqk(ie)%nbin=1
		end if
	end do
c	call chseqst(nor,tm3,ndr,eqk,ndeq,sta,ndst,mofeq,mofst)
	call datout1(nout,eqk,ndeq,sta,ndst)

	dt2=sqrt(dt2/ir)
	write(*,*)"Num of events/stations:",ndeq,ndst
	write(*,*)"Num of rays: total/chosed:",ia,ir
	write(*,*)"rms of residual without weighting",dt2
	end
c-------------------------------------------------------------------
	subroutine bingrp8azm(jaz,mbin0,mbin,minr,nbin,bin,nrbin)
	dimension bin(2,mbin),jaz(mbin0),nrbin(mbin)
	if(mod(mbin0,mbin).ne.0)then
	write(*,*)"The Mbin0 should be the times of Mbin in bingrp8azm"
	stop
	end if
	nbin=0
	wbin=360.0/mbin0
		nb0=mbin0/mbin		
		nr=0		
		j0=1		
		j1=1
		do i=1,mbin
			do j=(i-1)*nb0+1,i*nb0
				nr=nr+jaz(j)
				j1=j1+1
			end do
			if(nr.ge.minr)then
				nbin=nbin+1
				nrbin(nbin)=nr
				bin(1,nbin)=(j0-1)*wbin		
				bin(2,nbin)=(j1-1)*wbin
				if(nbin.eq.mbin.and.bin(2,nbin).ne.360.0)then
					bin(2,nbin)=360.0
					nrbin(nbin)=nr+jaz(mbin0)
					goto 886
				end if		
				nr=0; 		j0=j1
			elseif(i.eq.mbin.and.nr.gt.0)then
				if(nbin.eq.0)nbin=1
				nrbin(nbin)=nrbin(nbin)+nr
				bin(2,nbin)=360.0
			end if
886		end do
	end
c-------------------------------------------------------------------
	subroutine bingrp8rnoazm(jaz,mbin0,mbin,minr,nbin,bin,nrbin)
	dimension bin(2,mbin),jaz(mbin0),nrbin(mbin)
	wbin=360/mbin0
		nr=0; nbin=0
		do j=1,mbin0
			nr=nr+jaz(j)
		end do
		n1bin=nr/mbin
		if(n1bin.lt.minr)n1bin=minr
		n0=0; j0=1; j1=1
		do j=1,mbin0
			n0=n0+jaz(j)
			if(n0.ge.n1bin)then
				j1=j
				nbin=nbin+1
				nrbin(nbin)=n0
				bin(1,nbin)=(j0-1)*wbin		
				bin(2,nbin)=j1*wbin
				nrbin(nbin)=n0
				if(nbin.eq.mbin.and.bin(2,nbin).ne.360.0)then
					bin(2,nbin)=360.0
					nrbin(nbin)=n0+jaz(mbin0)
					return
				end if		
				n0=0; j0=j+1; j1=j+1
			end if
		end do
	end
c-------------------------------------------------------------------
	function vapp(edp,dmh,vcr,dst,tim)
		del=dst/180.0*3.1415926*(6371.0-dmh)
		tapp=tim-(dmh-edp)/10.0-dmh/10.0
		vapp=del/tapp
	end
c-------choose data by record number of every earthquake and station---------
	subroutine chseqst(nor,tm3,ndr,eqk,ndeq,sta,ndst,meq,mst)
	include "datatype.inc"
	dimension nor(2,ndr),tm3(3,ndr)
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)
	ndscd=0
222	nnn=0
	!---1--------cancel the events with less records than meq
	neq0=0
	do ie=1,ndeq
		if(eqk(ie)%nrd.lt.meq)then
			nr0=0
			do ir=1,ndr
				if(nor(1,ir).eq.ie)then
					nnn=nnn+1
					ist=nor(2,ir)
					sta(ist)%nrd=sta(ist)%nrd-1
				else
					nr0=nr0+1
					tm3(1,nr0)=tm3(1,ir);tm3(2,nr0)=tm3(2,ir)
					nor(1,nr0)=nor(1,ir)
					nor(2,nr0)=nor(2,ir)
				end if
			end do
			ndr=nr0
		else
			neq0=neq0+1
			do ir=1,ndr
				if(nor(1,ir).eq.ie)nor(1,ir)=neq0
			end do			
			eqk(neq0)=eqk(ie)
c			eqk(neq0).nrd=eqk(ie).nrd
c			eqk(neq0).elo=eqk(ie).elo
c			eqk(neq0).ela=eqk(ie).ela
c			eqk(neq0).edp=eqk(ie).edp
		end if
	end do
	ndeq=neq0
	!---2--------cancel the stations with less records than mst
	nst0=0
	do is=1,ndst
		if(sta(is)%nrd.lt.mst)then
			nr0=0
			do ir=1,ndr
				if(nor(2,ir).eq.is)then
					nnn=nnn+1
					ieq=nor(1,ir)
					eqk(ieq)%nrd=eqk(ieq)%nrd-1
				else
					nr0=nr0+1
					tm3(1,nr0)=tm3(1,ir);tm3(2,nr0)=tm3(2,ir)
					nor(1,nr0)=nor(1,ir)
					nor(2,nr0)=nor(2,ir)
				end if
			end do
			ndr=nr0
		else
			nst0=nst0+1
			do ir=1,ndr
				if(nor(2,ir).eq.is)nor(2,ir)=nst0
			end do			
			sta(nst0)=sta(is)
c			sta(nst0).nrd=sta(is).nrd
c			sta(nst0).slo=sta(is).slo
c			sta(nst0).sla=sta(is).sla
c			sta(nst0).shg=sta(is).shg
		end if
	end do
	ndst=nst0
	ndscd=ndscd+nnn
	if(nnn.gt.0)goto 222
	end

c-----Initilize the model------------------------------
	subroutine inimod(nx,ny,nbk,mdva,par0)
	real mdva(nbk,3),par0(3)

	do ip=2,3
		do iy=1,ny
			id=(iy-1)*nx
			do ix=1,nx
				mdva(id+ix,ip)=0.0
			end do
		end do
	end do
	do iy=1,ny
		id=(iy-1)*nx
		do ix=1,nx
			mdva(id+ix,1)=par0(3)
		end do
	end do
	end
c-----------------------
	subroutine statistic(hstgrm,ncl,fbgn,fend,fun,i3)
	integer hstgrm(3,ncl)
	if(fun.ge.fbgn.and.fun.le.fend)then
		df=(fend-fbgn)/(ncl-1)
		nf=int((fun-fbgn+0.5*df)/df)+1
		if(nf.gt.ncl)nf=ncl
		hstgrm(i3,nf)=hstgrm(i3,nf)+1
	end if
	end
c--------------------------------------------------------------
	function ttime(ipar,mdva,nx,ny,nbr,nobk,segbk,azmbk)
	dimension nobk(nbr),segbk(nbr),azmbk(nbr)
	real mdva(nx*ny,3)
	tt=0.0
	do ir=1,nbr
		nb=nobk(ir)
		vk=mdva(nb,1)
		ak=mdva(nb,2)
		bk=mdva(nb,3)
		seg=segbk(ir)
		az=azmbk(ir)
c		write(*,*)"vk/ak/bk/az=",vk,ak,bk,az
		if(ipar.eq.0)then
			tt=tt+seg*(1.0/vk+ak*cos(2*az)+bk*sin(2*az))
		elseif(ipar.eq.1)then
			tt=tt+seg*(1.0/vk)
		end if
	end do
	ttime=tt
	end 
!-------------------------
	function inarea(area,po,pa)
	dimension area(3,2)
	inarea=0
	if(po.ge.area(1,1).and.po.le.area(2,1))then
		if(pa.ge.area(1,2).and.pa.le.area(2,2))then
			inarea=1
		end if
	end if
	end 
c--------------------------------------------------------
	subroutine SpikeMod(fmck,luinv,mdva,cmode,regs)
	include "surfpath.par"
	character str(3)*20,fmck*80,ss*80
	real mdva(nbk,3),cmode(6),regs(4,20)
	!cmode(1/2/3/4)==mode/width 0f pattern/# of 0-perturbation/offset
	str(1)="/pn Velocity"
	str(2)="/Ak of Anisotropy"
	str(3)="/Bk of Anisotropy"

	f2s=(re-par0(1))/re
	s2f=re/(re-par0(1))
	v0=par0(3);		
	dv=0.3*s2f; 
	a0=0.0;		da=0.00; 					!da=0.0
	b0=0.0;		db=0.00
	mdva=0.0

	if(cmode(1).eq.7)then
		nrg=int(cmode(2))
		dv=cmode(3)
		!Average model
		do iy=1,ny
			id=(iy-1)*nx
			do ix=1,nx
				mdva(id+ix,1)=v0
				mdva(id+ix,2)=a0
				mdva(id+ix,3)=b0
			end do
		end do
		!Setup spikes
		do ir=1,nrg
			nx0=int((regs(1,ir)-xo0)/dx)+1
			nx1=int((regs(2,ir)-xo0)/dx)
			ny0=int((regs(3,ir)-ya0)/dy)+1
			ny1=int((regs(4,ir)-ya0)/dy)
			if (nx0.lt.1) nx0=1
			if (nx1.gt.nx) nx1=nx
			if (ny0.lt.1) ny0=1
			if (ny1.gt.ny) ny1=ny
			!Spike model
			do iy=ny0,ny1
				id=(iy-1)*nx
				do ix=nx0,nx1
					mdva(id+ix,1)=v0+dv
				end do
			end do
		end do
	endif
	
	write(luinv,*)"//Check Objective Model"
c	write(luinv,*)v0
	vsp=v0*f2s
	write(luinv,*)vsp
	write(luinv,*)xo0,xo1,ya0,ya1
	write(luinv,*)nx,ny

	write(luinv,*)str(1)
	do iy=1,ny
		id=(iy-1)*nx
		write(luinv,'(500f9.5,1x)')(mdva(id+ix,1)*f2s,ix=1,nx)
	end do
	do ip=2,3
		write(luinv,*)str(ip)
		do iy=1,ny
			id=(iy-1)*nx
			write(luinv,'(500f9.5,1x)')(mdva(id+ix,ip)/f2s,ix=1,nx)
		end do
  	end do

	open(221,file="tempmodel.dat")
	do iy=1,ny
	yy=y0+(iy-1)*dy+dy/2.0
	do ix=1,nx
		xx=x0+(ix-1)*dx+dx/2.0
		id=(iy-1)*nx+ix
		call anisotropy(mdva(id,2)/f2s,mdva(id,3)/f2s,amp,ang)
		write(221,'(50f11.5,1x)')xx,yy,mdva(id,1)*f2s, amp,ang
	end do
	end do
	close(221)
	end
c--------------------------------------------------------
	subroutine ckmod(fmck,luinv,mdva,cmode)
	include "surfpath.par"
	character str(3)*20,fmck*80,ss*80
	real mdva(nbk,3),cmode(6),regs(4,20)
	!cmode(1/2/3/4)==mode/width 0f pattern/# of 0-perturbation/offset
	str(1)="/pn Velocity"
	str(2)="/Ak of Anisotropy"
	str(3)="/Bk of Anisotropy"

	f2s=(re-par0(1))/re
	s2f=re/(re-par0(1))
	v0=par0(3);		
	dv=0.3*s2f; 
	a0=0.0;		da=0.00; 					!da=0.0
	b0=0.0;		db=0.00
	mdva=0.0
	if(cmode(1).eq.1)then  !boxcar mode
		ncx=(cmode(2)+cmode(3))/dx; 		dpx=3.1415926/2.0/ncx;
		ncy=(cmode(2)+cmode(3))/dy; 		dpy=3.1415926/2.0/ncy;
		sgny=1.0;
		da=0.0
		db=0.015
		do iy=1,ny
			sgnx=1.0; 		
			s0=1.0;
   			iyy=iy+cmode(5)-1
c			if(mod(iyy,int(ncy)).ge.cmode(2))s0=0.0
			if(mod(iyy,int(ncy)).eq.0)sgny=sgny*(-1.0)
			id=(iy-1)*nx
			do ix=1,nx
				s1=1.0
				ixx=ix+cmode(4)-1
c				if(mod(ixx,int(ncx)).ge.cmode(2))s1=0.0
				if(mod(ixx,int(ncx)).eq.0)sgnx=sgnx*(-1.0)
				!model of velocity
				mdva(id+ix,1)=v0+dv*sgnx*sgny*s0*s1*cmode(6)
				!model of anisotropy
				mdva(id+ix,2)=a0+da*sgnx*sgny*s0*s1*cmode(6)
c				mdva(id+ix,2)=a0+da
				mdva(id+ix,3)=b0+db*sgnx*sgny*s0*s1*cmode(6)
			end do
		end do
	elseif(cmode(1).eq.2.or.cmode(1).eq.3)then  !sinusoidal mode
		dv=0.3*s2f; 
		da=0.00625;  !2%
 		da=0.015;    !5%
		if (cmode(3).eq.2.0)then
			da=0.0
		elseif (cmode(3).eq.3.0)then
			dv=0.0
		end if
		wptn=cmode(2)
		do iy=1,ny
			id=(iy-1)*nx
			do ix=1,nx
				mdva(id+ix,1)=v0+dv*sinck(ix,iy,dx,dy,wptn)	
				!-----------checkerboard model for anisotropy
				mdva(id+ix,2)=a0+da*sinck(ix,iy,dx,dy,wptn)	
				mdva(id+ix,3)=b0+db*sinck(ix,iy,dx,dy,wptn)	
			end do
		end do
	elseif(cmode(1).eq.4.or.cmode(1).eq.5.or.cmode(1).eq.6)then !real model
	nreg=1
c	set the regions to be excluded
c	regs(1,1)=100.0;	regs(2,1)=105.0;  	
c	regs(3,1)=32.0;  	regs(4,1)=41.0;
c	regs(1,2)=		regs(2,2)=  	regs(3,2)=  	regs(4,2)=
c	regs(1,3)=		regs(2,3)=  	regs(3,3)=  	regs(4,3)=
c	regs(1,4)=		regs(2,4)=  	regs(3,4)=  	regs(4,4)=

		lum=221
		dv=0.3*s2f; da=0.00375/s2f; 		
		open(lum,file=fmck,status="old")
11		read(lum,"(50a)")ss
		if(ss(2:7).ne."//Real")goto 11
		read(lum,*)v0
		read(lum,*)x0,x1,y0,y1
		read(lum,*)nx0,ny0
		if(nx0.ne.nx.or.ny0.ne.ny)then
		   write(*,*)"Input Model's grid-size differ with the inversion!"
			stop	
		end if
		read(lum,*)
		do j=1,ny
			read(lum,*)(mdva((j-1)*nx+i,1),i=1,nx)
		end do
		if(cmode(1).eq.5)goto 887
101		read(lum,"(50a)",end=884)ss
		if(ss(2:4).ne."/Ak")goto 101
		goto 885
884		write(*,*)"Can not read anisotropy,Not included in inversion"
		stop
885		do j=1,ny
			read(lum,*)(mdva((j-1)*nx+i,2),i=1,nx)
		end do
		read(lum,*)
		do j=1,ny
			read(lum,*)(mdva((j-1)*nx+i,3),i=1,nx)
		end do
887		close(lum)
		do iy=1,ny
			id=(iy-1)*nx
			do ix=1,nx
				mdva(id+ix,1)=mdva(id+ix,1)*s2f
				mdva(id+ix,2)=mdva(id+ix,2)/s2f
				mdva(id+ix,3)=mdva(id+ix,3)/s2f
			end do
		end do

		v00=v0*s2f
		do iy=1,ny
			id=(iy-1)*nx
			do ix=1,nx
				xlo=x0+(ix-1)*dx
				yla=y0+(iy-1)*dy
c				kin=inregs(regs,nreg,xlo,yla)
				kin=0
				v01=v00-cmode(2)-0.05
				v02=v00+cmode(2)
				if(kin.eq.1)then
c					mdva(id+ix,1)=v00
				elseif(mdva(id+ix,1).lt.v01)then
					mdva(id+ix,1)=v00-dv
				elseif(mdva(id+ix,1).gt.v02)then
					mdva(id+ix,1)=(v00+dv)
				else
					mdva(id+ix,1)=v00
				end if
				!-----------checkerboard model for anisotropy
				call anisotropy(mdva(id+ix,2),mdva(id+ix,3),amp,ang)
				if(amp.lt.cmode(3).or.cmode(1).eq.5)then
					mdva(id+ix,2)=0.0
					mdva(id+ix,3)=0.0
				elseif(amp.gt.cmode(3))then
					aaa=da
					!aaa=amp*0.125
					call inverse_anisotropy(aaa,-ang,ak,bk)
					!write(*,*)mdva(id+ix,2),ak,mdva(id+ix,3),bk
c     					mdva(id+ix,2)=ak
c					mdva(id+ix,3)=bk
				end if
				if(cmode(1).eq.6)mdva(id+ix,1)=v0
			end do
		end do
	else
	endif
	
	write(luinv,*)"//Check Objective Model"
c	write(luinv,*)v0
	vsp=v0*f2s
	write(luinv,*)vsp
	write(luinv,*)xo0,xo1,ya0,ya1
	write(luinv,*)nx,ny

	write(luinv,*)str(1)
	do iy=1,ny
		id=(iy-1)*nx
		write(luinv,'(500f9.5,1x)')(mdva(id+ix,1)*f2s,ix=1,nx)
	end do
	do ip=2,3
		write(luinv,*)str(ip)
		do iy=1,ny
			id=(iy-1)*nx
			write(luinv,'(500f9.5,1x)')(mdva(id+ix,ip)/f2s,ix=1,nx)
		end do
  	end do

	open(221,file="tempmodel.dat")
	do iy=1,ny
	yy=y0+(iy-1)*dy+dy/2.0
	do ix=1,nx
		xx=x0+(ix-1)*dx+dx/2.0
		id=(iy-1)*nx+ix
		call anisotropy(mdva(id,2)/f2s,mdva(id,3)/f2s,amp,ang)
		write(221,'(50f11.5,1x)')xx,yy,mdva(id,1)*f2s, amp,ang
	end do
	end do
	close(221)
	end
c--------------------------------------------------------------------
	function inregs(regs,nreg,xlo,yla)
	dimension regs(4,nreg)
	inregs=0
	do ireg=1,nreg
		if(xlo.ge.regs(1,ireg).and.xlo.le.regs(2,ireg).and.
     +				yla.ge.regs(3,ireg).and.yla.le.regs(4,ireg))then
			inregs=1
			return		
		end if
	end do
	end
c--------------------------------------------------------------------
	subroutine anisotropy(aij,bij,amp,ang)
	pp=180.0/3.1415926
		am=sqrt(aij*aij+bij*bij)
		!------------direction------------------
			if(abs(aij).le.0.0000001.and.abs(aij).le.0.0000001)then
				an=0.0; am=0.0
			elseif(abs(aij).le.0.0000001)then
				an=90.0
				if(bij.gt.0.0)an=270.0
c			elseif(abs(bij).le.0.000001)then
c				an=0.0;
c				if(aij.gt.0.0)an=180.0
			else
				an=atan(bij/aij)*pp
				if(aij.lt.0.0.and.bij.lt.0.0)then
					an=an+0.0
				elseif(aij.gt.0.0.and.bij.lt.0.0)then
					an=an+180.0
				elseif(aij.gt.0.0.and.bij.gt.0.0)then
					an=an+180.0
				elseif(aij.lt.0.0.and.bij.gt.0.0)then
					an=an+360.0
				end if
			end if


		ang=-an/2.0		  !!negative for plotting purpose?
		amp=am/0.125
	end
c--------------------------------------------------------------------
	subroutine inverse_anisotropy(amp,ang,ak,bk)
		pp=180.0/3.1415926
		if(ang.lt.0.0)ang=ang+180.0
		rang=ang/pp
		if(ang.ne.45.0.and.ang.ne.135.0)then
			aij=amp/sqrt(1.0+tan(2.0*rang)*tan(2.0*rang))
			bij=sqrt(amp*amp-aij*aij)
		end if
		!------------direction------------------
			if(amp.le.0.0000001)then
				ak=0.0; bk=0.0
			elseif(ang.eq.45.0)then
				bk=-amp; ak=0.0
			elseif(ang.eq.135.0)then
				bk=amp; ak=0.0
			elseif(ang.ge.0.0.and.ang.lt.45.0)then
				ak=-aij; bk=-bij
			elseif(ang.gt.45.0.and.ang.lt.90.0)then
				ak=aij; bk=-bij
			elseif(ang.ge.90.0.and.ang.lt.135.0)then
				ak=aij; bk=bij
			elseif(ang.gt.135.0.and.ang.lt.180.0)then
				ak=-aij; bk=bij
			end if
	end
!------------calculation of sinusoidal checkerboard------
	function sinck(ix,iy,dx,dy,dlen)
	pi=3.1415926
	ddx=dlen/dx;
	dd=pi/ddx
	xx=ix*dd-dd/2.0
	shift=0.0
	if(mod(ddx,real(int(ddx))).eq.0.0)shift=pi/2.0-dd/2.0
	xx=xx+shift

	ddy=dlen/dy;
	dd=pi/ddy
	yy=iy*dd-dd/2.0
	shift=0.0
	if(mod(ddy,real(int(ddy))).eq.0.0)shift=pi/2.0-dd/2.0
	yy=yy+shift

	sinck=sin(xx)*sin(yy)
	end
C----------------------------------------------
	subroutine smooth(vv,nx,ny,fsmth)
	dimension vv(nx,ny),vo(nx,ny)
	do 100 ix=1,nx
	do 100 iy=1,ny
		nn=0
		va=0.0
		ww=0.0
		do 101 i0=ix-5,ix+5
		do 101 i1=iy-5,iy+5
			if(i0.ge.1.and.i0.le.nx.and.i1.ge.1.and.i1.le.ny)then
				if(i0.eq.ix.and.i1.eq.iy)goto 101
				nn=nn+1
				x0=i0		; y0=i1
				x1=ix		; y0=iy
				wi=weight(x0,y0,x1,y1)
				va=va+vv(i0,i1)*wi
				ww=ww+wi
			end if
101		end do
		dv=va/ww-vv(ix,iy)
		dv=dv*fsmth
		vv(ix,iy)=vv(ix,iy)+dv
100		end do

		do ix=1,nx
			do iy=1,ny
				vv(ix,iy)=vv(ix,iy)
			end do
		end do
		end

		function weight(x0,y0,x1,y1)
			xx=(x0-x1)*(x0-x1)
			yy=(y0-y1)*(y0-y1)
			dis=sqrt(xx+yy)
			weight=1.0/dis
	end
