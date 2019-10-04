	program surfinv
c	use portlib
c	use	msflib
	include "surfpath.par"
	include "datatype.inc"
	data area/xo0,xo1,dx,ya0,ya1,dy/
	character pth*80,frc*80,finv*80,fdat*80,fmck*80,head*80,str*6
	dimension tim(3,nry),nor(2,nry),qual(nry)
	dimension norb(2,nry),timb(3,nry),qualb(nry),idxb(nry)
	integer(4) irandom
	type (station) sta(nst)
	type (earthquake) eqk(neq)
	type (partom) ipar

	open(1,file="surftom.ini",status="old")
	read(1,*)head; write(*,*)head
	read(1,*)pth; write(*,*)pth
	read(1,*)frc; write(*,*)frc
	read(1,*)(ipar%cmode(i),i=1,3) 
	read(1,*)nboot
	read(1,*)ipar%isplit,ipar%nr1gp
	read(1,*)ipar%mbin0,ipar%mbst,ipar%minrst, ipar%mbeq, ipar%minreq
	read(1,*)ipar%icrr
        read(1,*)ipar%ipar  !0--velo+aniso+st+eq//1---no anisotropy
	read(1,*)ipar%idv  !0--ds inversion  1--v inversion 2--dv
	read(1,*)ipar%iwgh,ipar%iwdat  !1-No weighting
	read(1,*)ipar%iterr,ipar%terr,ipar%terrdev !add random error to time
	read(1,*)ipar%vsmth, ipar%asmth  !smoothing factors//asmth
	read(1,*)(wgh(i),i=1,4)  !//weights for Velocity/anis/stations/earthquake
	read(1,*)kdt,dtmin,dtmax  !residual limit
	read(1,*)kds,dsmin,dsmax  !distance limit
	read(1,*)ipar%itmax  !Maximum iteration times
	read(1,*)finv
	read(1,*)fdat
	if (ipar%cmode(1).ge.4.and.ipar%cmode(1).le.6)then
		read(1,*)
		read(1,*)fmck
		fmck=trim(pth)//trim(fmck)
	end if
	if (ipar%cmode(1).eq.7)then
11		read (1,"(a6)")str
		if(str.ne."/Spike")goto 11
		read(1,*)ipar%cmode(2),ipar%cmode(3)
		do i1=1,int(ipar%cmode(2)) 
			read(1,*)(ipar%subreg(i2,i1),i2=1,4)
		enddo
	end if
	close(1)
	write(*,*)"fdat=",fdat
	!---------------------Data Preparing-----------------
	fdat=trim(pth)//trim(fdat)
	finv=trim(pth)//trim(finv)
	write(*,*)"fdat=",fdat,"finv=",finv
	lulog=90	
	open(lulog,file="invlog.dat")
	luout=111  
	open(luout,file=fdat)
	luinv=112	
	open(luinv,file=finv)
c	---Output inversion parameters----------
      call datout0(luout,head,ipar,frc,finv,fdat)
      call datout0(luinv,head,ipar,frc,finv,fdat)
      write(*,*)"Now start computing--"
	!	--Real inversion--
	if(ipar%cmode(1).eq.0)then
	CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     + ndr,eqk,sta,ipar)
		write(*,*)"//Real Inversion--Output Model of Inversion"
		write(luinv,*)"//Real Inversion"
		call inversion(luout,luinv,lulog,
     +			nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
		do ib=1,nboot
			write(*,*)"//iboot=",ib,"/",nboot
			write(luinv,*)"//iboot=",ib,"/",nboot
			!generate the random index number of rays
			ndrb=int(1.0*ndr);
			do ir=1,ndrb
101				irandom=irand(0)
				idxr=mod(irandom,ndrb)+1
				do jb=1,ir-1
					!if (idxr.eq.idxb(jb))goto 101
				end do
				norb(1,ir)=nor(1,idxr)
				norb(2,ir)=nor(2,idxr)
				qualb(ir)=qual(idxr)
				timb(1,ir)=tim(1,idxr)
				timb(2,ir)=tim(2,idxr)
				idxb(ir)=idxr;
				if(mod(idxr,1000).eq.0)then
					write(lulog,*)idxr,timb(1,ir),timb(2,ir)
				end if
			end do
			call inversion(luout,luinv,lulog,
     +			norb,qualb,eqk,sta,ndeq,ndst,ndrb,timb,ipar,1)
		end do
	!	--Real inversion & Checkerboard--
      elseif(ipar%cmode(1).eq.1)then
c	----Real inversion first
		ipar%cmode(1)=0		
		CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     +		ndr,eqk,sta,ipar)
		write(*,*)"//Real Inversion--Output Model of Inversion"
		write(luinv,*)"//Real Inversion"
		call inversion(luout,luinv,lulog,
     +			nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
c	---checkerboard test--
c	---nck0=# of checkerboard test along one direction. total # of tests==nck0*nck0
c	---
		ipar%cmode(1)=1		
		nck0=1
c		nck0=int((ipar%cmode(2)+ipar%cmode(3))/ipar%cmode(2)+0.5)
		ipar%cmode(6)=1.0
		do 111 ick=1,nck0
		ipar%cmode(4)=(ick-1)*ipar%cmode(2)
		do 111 jck=1,nck0
		ipar%cmode(5)=(jck-1)*ipar%cmode(2)
		ipar%cmode(6)=ipar%cmode(6)*(-1.0)
		nck=(ick-1)*nck0+jck
		write(*,*)"//Inversion  ",nck,"/",nck0*nck0
		CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     +		ndr,eqk,sta,ipar)
		write(*,*)"//Check Inversion--Output Model of Inversion"
		write(luinv,*)"//Check Inversion  ",nck,"/",nck0*nck0
		call inversion(luout,luinv,lulog,
     +		nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
111		enddo
	!Sinusoidal Boxcar checkerboard
	elseif(ipar%cmode(1).eq.2)then
		CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     +		ndr,eqk,sta,ipar)
		write(*,*)"//Check Inversion--Output Model of Inversion"
		write(luinv,*)"//Check Inversion  "
		call inversion(luout,luinv,lulog,
     +		nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
	elseif(ipar%cmode(1).eq.5)then
		CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     +		ndr,eqk,sta,ipar)
		write(*,*)"//Check Real Inversion--Output Model of Inversion"
		write(luinv,*)"//Real Inversion  "
		call inversion(luout,luinv,lulog,
     +		nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
	elseif(ipar%cmode(1).eq.7)then
		CALL surfdata(luinv,fmck,frc,luout,nor,tim,qual,ndeq,ndst,
     +		ndr,eqk,sta,ipar)
		write(*,*)"//Spike Check"
		write(luinv,*)"//Real Inversion  "
		call inversion(luout,luinv,lulog,
     +		nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,0)
	end if

	close(luinv)
	close(lulog)
	end


c-------------------random gaussian errors...
	subroutine random_gaussian(scut,dev,srand)
100		CALL RANDOM_NUMBER (ran1)
		CALL RANDOM_NUMBER (ran2)
		ran1=2.0*ran1-1.0
		ran2=2.0*ran2-1.0
		r=ran1*ran1+ran2*ran2
		if (r.ge.1.0.or.r.eq.0.0) goto 100		
		fac=sqrt(-2.0*log(r)/r)*dev  !*dev---added myself, not proved.
		r1=ran1*fac
		r2=ran2*fac

		if(abs(r1).le.scut)then
			srand=r1
		elseif(abs(r2).le.scut)then
			srand=r2
		else
			goto 100
		end if
	end
