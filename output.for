c-------Output inversion parameters-----------------------
	subroutine datout0(lu,head,ipar,fdat,finv,foth)
	include "surfpath.par"
	include "datatype.inc"
	type (partom)::ipar
	character*80 head,fdat,finv,foth
	write(lu,'(a4 a80)')" ///",head
	write(lu,'(a11 a80)')" Data File:",fdat
	write(lu,'(a19 a80)')" Output Model File:",finv
	write(lu,'(a19 a80)')" Output Info. File:",foth
	write(lu,*)"Area:Lon1 Lon2 DLon Lat1 Lat2 DLat"
	write(lu,'(6f8.3)')xo0,xo1,dx,ya0,ya1,dy
	write(lu,*)"Inversion Parameters:"
	write(lu,*)"Time (Residual) Cut:"
	write(lu,'(I3,2f8.3)')kdt,dtmin,dtmax
	write(lu,*)"Smoothing (Vel & Ans):"
	write(lu,'(2f8.3)')ipar%vsmth, ipar%asmth
	write(lu,*)"Data weighting:"
	write(lu,*)ipar%iwgh,ipar%iwdat
	write(lu,'(4f8.3)')(wgh(i),i=1,4) 
	write(lu,*)"Distance Cut (deg):"
	write(lu,'(I3,2f8.3)')kds,dsmin,dsmax
	if (ipar%cmode(1).eq.7)then
		write(lu,*)"/Spike Check: Locations of spikes:"
		write(lu,*)ipar%cmode(2),ipar%cmode(3)
		do i1=1,int(ipar%cmode(2)) 
			write(lu,*)(ipar%subreg(i2,i1),i2=1,4)
		enddo
	end if
	write(lu,*)"///End of parameter list///" 
	write(lu,*)
	return

c	read(1,*)(ipar.cmode(i),i=1,3) !0--real inversion//999--checkerboard 
c	read(1,*)nboot
c	read(1,*)ipar.icrr
c c     read(1,*)ipar.ipar  !0--velo+aniso+st+eq//1---no anisotropy
c	read(1,*)ipar.idv  !0--ds inversion  1--v inversion 2--dv
c	read(1,*)ipar.iterr,ipar.terr,ipar.terrdev !add random error to time
c	read(1,*)ipar.itmax  !Maximum iteration times

	end
c-------Output events and station information-----------------------
	subroutine datout1(nout,eqk,ndeq,sta,ndst)
	include "datatype.inc"
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)

	write(nout,*)"/events used for the pn-tomography"
	write(nout,*)ndeq
	do i=1,ndeq
		write(nout,"(f8.3,1x,f7.3,1x,f6.1,1x,i4,1x,a5)")
     +	eqk(i)%elo,eqk(i)%ela,eqk(i)%edp,eqk(i)%nrd,eqk(i)%eqn
	end do

	write(nout,*)"/stations used for the pn-tomography"
	write(nout,*)ndst
	do i=1,ndst
		write(nout,"(f8.3,1x,f7.3,1x,f6.1,1x,i4,1x,a5)")
     +	sta(i)%slo,sta(i)%sla,sta(i)%shg,sta(i)%nrd,sta(i)%stn
	end do
	return
	end
c-------output events and station information----------------------------
	subroutine datout2(nout,nx,ny,rcov)
	dimension rcov(4,nx*ny)

	write(nout,*)"/Numbers of rays In every block"
	write(nout,*)nx,ny
	do iy=1,ny
		write(nout,"(100(i5,1x))")(int(rcov(4,(iy-1)*nx+ix)),ix=1,nx)
	end do

	write(nout,*)"/Length of rays In every block"
	write(nout,*)nx,ny
	do iy=1,ny
		write(nout,99)(rcov(1,(iy-1)*nx+ix)/111.199,ix=1,nx)
	end do

	write(nout,*)"/EW-Length of rays in every block"
	write(nout,*)nx,ny
	do iy=1,ny
		write(nout,99)(rcov(2,(iy-1)*nx+ix)/111.199,ix=1,nx)
	end do

	write(nout,*)"/NS-Length of every block in N-S Direction"
	write(nout,*)nx,ny
	do iy=1,ny
		write(nout,99)(rcov(3,(iy-1)*nx+ix)/111.199,ix=1,nx)
	end do
99	format(200(f8.3))
	end
C-----------------------------------------
	subroutine datout3(nout,ndr,nor,tm3,qual)
	parameter(nhgrm=41,	dtdw=-10, dtup=10)
	dimension nor(2,ndr),tm3(3,ndr),qual(ndr)
	integer hstgrm(3,nhgrm)
	write(nout,*)"/rays used for the pn-tomography"
		write(nout,*)ndr
		do i=1,ndr
			dt12=(tm3(1,i)-tm3(2,i))*qual(i)
			dt13=(tm3(1,i)-tm3(3,i))*qual(i)
			dt23=(tm3(3,i)-tm3(2,i))*qual(i)
			write(nout,90)nor(1,i),nor(2,i),tm3(1,i),tm3(2,i),tm3(3,i)
     +				,dt12,dt13,dt23
			call statistic(hstgrm,nhgrm,dtdw,dtup,dt12,1)
			call statistic(hstgrm,nhgrm,dtdw,dtup,dt13,2)
			call statistic(hstgrm,nhgrm,dtdw,dtup,dt23,3)
			avg12=avg12+dt12
			avg13=avg13+dt13
			dev12=dev12+dt12*dt12
			dev13=dev13+dt13*dt13
		end do
90	format(2(1x,i5),3(1x,f8.3),3(1x,f8.4))

	avg12=avg12/ndr
	avg13=avg13/ndr
	rms12=sqrt(dev12/ndr)
	rms13=sqrt(dev13/ndr)
	dev12=sqrt(dev12/ndr-avg12*avg12)
	dev13=sqrt(dev13/ndr-avg13*avg13)

	write(nout,*)"/Residual Distribution statistic"
	write(nout,*)"rms/avg/dev before inversion", rms12,avg12,dev12
	write(nout,*)"rms/avg/dev after  inversion", rms13,avg13,dev13
	write(*,*)"rms/avg/dev before inversion(w)", rms12,avg12,dev12
	write(*,*)"rms/avg/dev after  inversion", rms13,avg13,dev13
	write(nout,*)nhgrm
	dt=(dtup-dtdw)/(nhgrm-1)
	do ih=1,nhgrm
		tt=(ih-1)*dt+dtdw
		write(nout,*)tt,hstgrm(1,ih),hstgrm(2,ih),hstgrm(3,ih)
	end do
	end
!----------------------out put--the inversion result-------------
	subroutine invout(luinv,resu,ipar,sta,ndst,eqk,ndeq,s0,a0,b0)
	include "surfpath.par"
	include "datatype.inc"
	dimension resu(npar),dt(100)
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)
	type (partom)::ipar

	if(ipar%ipar.eq.0)then  !include anisotropy
		np=3
	elseif(ipar%ipar.eq.1)then !don't include anisotropy
		np=1
	end if
c	ncol=np*nbk+ndst+ndeq   !number of colums of the matrix
	ncol=np*nbk   !number of colums:No st/eq params for surf wave inversion
	v0=1.0/s0
	!c real inversion	!output the inversion model
	f2s=(re-par0(1))/re
	vs0=v0*f2s
	ds0=re-re*exp(-par0(1)/re)

	write(luinv,*)vs0
	write(luinv,*)xo0,xo1,ya0,ya1
	write(luinv,*)nx,ny
	write(luinv,*)"/Pn Variation"
	do ia=1,ny
		il=(ia-1)*nx*np
		if(ipar%idv.eq.0)then
			do io=1,nx
				if(abs(resu(il+(io-1)*np+1)+s0).le.0.05)then
					resu(il+(io-1)*np+1)=0.0
				end if
			end do
		   write(luinv,"(200f9.4)")(f2s/
     +				   (resu(il+(io-1)*np+1)+s0),io=1,nx)
		elseif(ipar%idv.eq.1)then
			do io=1,nx
				if(abs(resu(il+(io-1)*np+1)).le.0.1)then
					resu(il+(io-1)*np+1)=s0
				end if
			end do
		  write(luinv,"(200f9.4)")(f2s/(resu(il+(io-1)*np+1)),io=1,nx) 	
		elseif(ipar%idv.eq.2)then
		   write(luinv,"(200f9.4)")(f2s*
     +				   (resu(il+(io-1)*np+1)+v0),io=1,nx)
		end if
	end do

	!output the anisotropy
	if(ipar%ipar.eq.0)then
		write(luinv,*)"/Ak of anisotropy"
		do ia=1,ny
		il=(ia-1)*nx*3
          write(luinv,"(200f9.5)")
     +		((resu(il+(io-1)*np+2)+a0)/f2s,io=1,nx)
		end do
		write(luinv,*)"/Bk of anisotropy"
		do ia=1,ny
		il=(ia-1)*nx*3
		write(luinv,"(200f9.5)")
     +		((resu(il+(io-1)*np+3)+b0)/f2s,io=1,nx)
		end do
	end if

887	continue
224	format(f8.3,1x,f7.3,1x,f7.3,1x,i4,f5.2,1x,f7.3,1x,i1,1x,i1)
	end
!----------------------out put--the inversion result-------------
	subroutine invoutold(luinv,resu,terr,dterr,eerr,derr,deerr,
     +				ipar,ityp,sta,ndst,eqk,ndeq,s0,a0,b0)
	include "surfpath.par"
	include "datatype.inc"
	dimension resu(npar),dterr(npar),deerr(npar)
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)

	if(ipar.eq.0)then  !include anisotropy
		np=3
	elseif(ipar.eq.1)then !don't include anisotropy
		np=1
	end if
	ncol=np*nbk+ndst+ndeq   !number of colums of the matrix
	v0=1.0/s0
	!c real inversion	!output the inversion model
	f2s=(re-par0(1))/re
	write(*,*)"s0=",s0
	vs0=v0*f2s
	ds0=re-re*exp(-par0(1)/re)
	write(*,*)"f2s=",f2s
	write(*,*)"vs0=",vs0," ds0=",ds0

	write(luinv,*)vs0,ds0
	write(luinv,*)xo0,xo1,ya0,ya1
	write(luinv,*)nx,ny
	write(luinv,*)"/Pn Variation"
	do ia=1,ny
		il=(ia-1)*nx*np
		if(ityp.eq.0)then
		   write(luinv,"(200f9.4)")(f2s/
     +				   (resu(il+(io-1)*np+1)+s0),io=1,nx)
		else
			do io=1,nx
				if(abs(resu(il+(io-1)*np+1)).le.0.1)then
					resu(il+(io-1)*np+1)=s0
				end if
			end do
		  write(luinv,"(200f9.4)")(f2s/(resu(il+(io-1)*np+1)),io=1,nx) 	
		end if
	end do
	!output the anisotropy
	if(ipar.eq.0)then
		write(luinv,*)"/Ak of anisotropy"
		do ia=1,ny
		il=(ia-1)*nx*3
          write(luinv,"(200f9.5)")
     +		((resu(il+(io-1)*np+2)+a0)/f2s,io=1,nx)
		end do
		write(luinv,*)"/Bk of anisotropy"
		do ia=1,ny
		il=(ia-1)*nx*3
		write(luinv,"(200f9.5)")
     +		((resu(il+(io-1)*np+3)+b0)/f2s,io=1,nx)
		end do
	end if
c t-check--!output the random error for check of traverl-time differenvce
	if(terr.eq.0.0)goto 886
	dmax=-111.0
	do ia=1,ny
		il=(ia-1)*nx*np
		do io=1,nx
			verr=f2s/(dterr(il+(io-1)*np+1)+s0)
			vinv=f2s/(resu(il+(io-1)*np+1)+s0)
			dterr(il+(io-1)*np+1)=verr-vinv
			de0=verr-vinv
			if(abs(de0).ge.10.0)then
				de0=0.0
				dterr(il+(io-1)*np+1)=de0
				write(*,*)io,ia,verr,vinv-verr
			end if
			if(dmax.lt.abs(de0))dmax=de0
			de=de+de0*de0
		end do
	end do
	de=sqrt(de/(nx*ny-2))
	write(luinv,*)"/Vt-err-Random error test:rmt=",de,dmax
	do ia=1,ny
		il=(ia-1)*nx*3
		write(luinv,"(200f9.5)")((dterr(il+(io-1)*np+1)),io=1,nx)
	end do
c e-check--!output the random error check of location differenvce
886	dmax=-111.0
	if((derr+eerr).eq.0.0)goto 887
	de=0
	do ia=1,ny
		il=(ia-1)*nx*np
		do io=1,nx
			verr=f2s/(deerr(il+(io-1)*np+1)+s0)
			vinv=f2s/(resu(il+(io-1)*np+1)+s0)
			deerr(il+(io-1)*np+1)=verr-vinv
			de0=verr-vinv
			if(abs(de0).ge.10.0)then
				de0=0.0
				deerr(il+(io-1)*np+1)=de0
				write(*,*)io,ia,verr,vinv-verr
			end if
			if(dmax.lt.abs(de0))dmax=de0
			de=de+de0*de0
		end do
	end do
	de=sqrt(de/(nx*ny-2))
	write(luinv,*)"/Ve-err-Random error test:rmt=",de,dmax
	do ia=1,ny
		il=(ia-1)*nx*3
		write(luinv,"(200f9.5)")((deerr(il+(io-1)*np+1)),io=1,nx)
	end do

887	continue
223	format(f8.3,1x,f7.3,1x,f6.1,1x,f7.3,1x,i4,f5.2,1x,f6.3,1x,a5)
224	format(f8.3,1x,f7.3,1x,f7.3,1x,i4,f5.2,1x,f7.3,1x,i1,1x,i1)
225	format(f8.3,1x,f7.3,1x,f8.3,1x,f7.3,1x,i4,1x,f5.2,1x,f6.3,1x,2i2)
	end
