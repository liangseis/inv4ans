!-------------------------------------------------------------------!
!-----------invert the dataset and output the result----------------!
!-------------------------------------------------------------------!
	subroutine inversion(luout,luinv,lulog,
     +			nor,qual,eqk,sta,ndeq,ndst,ndr,tim,ipar,krun)
	include "surfpath.par"
	include "datatype.inc"
	dimension tim(3,nry),tinv(nry),nor(2,nry),qual(nry)
	dimension resui(npar),resu(npar),dterr(npar),deerr(npar)
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)
	type (partom)::ipar

	npst=0
	do i=1,ndst
		npst=npst+sta(i)%nbin	
	end do
	npeq=0
	do i=1,ndeq
		npeq=npeq+eqk(i)%nbin	
	end do

c2	!-------step 2---------Form Jacobe Matrix----A/Ax=b-------
	cmode=ipar%cmode(1);
	kpar=ipar%ipar;
	ityp=ipar%idv;   iwgh=ipar%iwgh; itmax=ipar%itmax
	terr=ipar%terr;		terdev=ipar%terrdev; iterr=ipar%iterr
	dperr=ipar%derr;    eperr=ipar%deerr	!eperr=ipar.deerr
	vsmth=ipar%vsmth;		asmth=ipar%asmth
	ieerr=0

118	write(*,*)"form jacobe matrix-----------"
	if (krun.ne.2)then
		call 
     +	jcbm(luout,ipar,nor,qual,ndr,eqk,ndeq,sta,ndst,nrow,nz,tim)
		bkw0=bkw;nrow0=nrow; nz0=nz; nzvl0=nzvl;
	end if
	!!combine the parameter weighting and length weighting
	nrow=nrow0
	nz=nz0
	nzvl=nzvl0
      bkw=bkw0
	call formwgh(nbk,npar,sta,ndst,eqk,ndeq,bkw,wgh,ipar)  
	if(iwgh.eq.1)then  !iwgh=1------no weights are applied!!
		do ip=1,npar
			bkw(ip)=1.0
		end do
	end if
	!Weighting the inversion matrix--------
	do inz=1,nz
		fact=1.0
		if(ipar%idv.eq.2)then  !idv=1------directly invert for dv!!
			v00=modva(icol,1)   !only good for inversion without anisotropy
			fact=-1.0/v00/v00  !has to modify for anisotropy
		end if
		icol=nzcl(inz)
		nzvl(inz)=nzvl(inz)*bkw(icol)*fact;
	end do
	!Weighting data--------
c	if (ipar.iwdat.eq.0)then
c		do ir=1,nrow
c			wdat=qual(ir);
c			tim(1,ir)=tim(1,ir)*wdat;
c			tim(2,ir)=tim(2,ir)*wdat;
c		end do
c	end	if
	!The degree of smoothing-too large will decrese resolution----
	ndata=nrow
	if(cmode.ne.999)then
		call rglrz(ipar%ipar,vsmth,asmth,nrow,nz)  
	endif
	!no regularization for checkboard
	if(nrow.gt.nry)stop "too many rows in the inversion matrix"
	if(nz.gt.nnz)stop "too many non-zero in the inversion matrix"
	damp=0.0
	drmin=100.0
	if(ipar%ipar.eq.0)then  !include anisotropy
		np=3
	elseif(ipar%ipar.eq.1)then !don't include anisotropy
		np=1
	end if
c	ncol=np*nbk+npst+npeq   !number of colums of the matrix
	ncol=np*nbk			   !number of colums of the matrix
	write(*,*)"Solve equations-----------"
	!------no random error invrsion-------------
	write(*,*)"No random errors inversion"
c3	!-------step 3---------Form the right side---b of Ax=b------
	kwd=ipar%iwdat
	if(iterr.eq.0) then
		call rhs(ityp,terr,terdev,ndata,nrow,qual,kwd,tim,tinv,nry)
	else
		call rhs(ityp,0.0,0.0,ndata,nrow,qual,kwd,tim,tinv,nry)
	end if
c4	!-------step 4---------Solve the Linear Equation system//Ax=b---
	write(*,*)"-------1-------"
	if(ieerr.eq.0)then
111		resu=0.0
		call LSQR(lulog,itmax,DAMP,nrow,ncol,tinv,drmin,resu,
     +                  nnzr,nzcl,nzvl,nry,nnz)
	else
		call LSQR(lulog,itmax,DAMP,nrow,ncol,tinv,drmin,deerr
     +                  nnzr,nzcl,nzvl,nry,nnz)
	end if

	write(*,*)"-------2-------"

	if(iterr.ne.0 .and. terr.ne.0.0)then   !do another inv with rand error
		call rhs(ityp,terr,terrdev,ndata,nrow,qual,kwd,tim,tinv,nry)
		iterr=iterr+1
		write(*,*)"-----random errors inversion-----",iterr
		call LSQR(lulog,itmax,DAMP,nrow,ncol,tinv,drmin,dterr
     +                  nnzr,nzcl,nzvl,nry,nnz)
	end if

c	!-------step 5---------Output the result------------------------
	s0=1.0/par0(3)
	call weightout(bkw,resu,ncol)
	call invout(luinv,resu,ipar,sta,ndst,eqk,ndeq,s0,a0,b0)
	
	write(*,*)"iterr=",iterr,"krun=",krun
	if(terr.ne.0.0) then
		call weightout(bkw,dterr,ncol)
		call invout(luinv,dterr,ipar,sta,ndst,eqk,ndeq,s0,a0,b0)
	endif
	if(krun.eq.0)then
		call rsdchg(ndr,tim,nor,eqk,ndeq,sta,ndst,resu,ityp,ipar)
		call datout3(luout,ndata,nor,tim,qual)
	endif
	end 
!-----------------right side of inversion equation system--------
	subroutine rhs(ityp,terr,terdev,ndat,nrow,qual,kwd,tim,tinv,nry)
	dimension tim(3,nry),tinv(nry),qual(nry)
	wdat=1.0
	if(ityp.eq.0.or.ityp.eq.2)then
		do ir=1,ndat
			if (kwd.eq.0)then
				wdat=qual(ir)
			endif
			tinv(ir)=(tim(1,ir)-tim(2,ir))*wdat  
			!real travel time - calculated time
		end do
		do i=ndat+1,nrow
			tinv(i)=0.0
		end do
	elseif(ityp.eq.1)then
		do ir=1,ndat
			if (kwd.eq.0)then
				wdat=qual(ir)
			endif
			tinv(i)=tim(1,ir)*wdat
		end do
		do i=ndat+1,nrow
			tinv(i)=0.0
		end do
	end if
	!calculate the standard deviation of the residul before inversion.
	aa=0.0
	av=0.0
	do i=1,ndat
		av=av+tinv(i) 
		aa=aa+ tinv(i)*tinv(i)
	end do
	av=av/ndat
	aa=aa/ndat
	a0=sqrt(aa-av*av)
	!add random errors to data
	write(*,*)"ter0,terdev=",terr,terdev
	if(terr.ne.0.0)then
		do i=1,ndat
			if(terdev.ne.0.0)then
				call random_gaussian(terr,terdev,rand)
			else
				CALL RANDOM_NUMBER (ranr)
				CALL RANDOM_NUMBER (rans)
				ranr=ranr*terr				    !random residual:<<terr
				if(rans.ge.0.5)rand=ranr*(-1.0) !the sign of the random error
			end if
			tinv(i)=tinv(i)+rand
		end do
	end if
	!calculate the standard deviation of the residul before inversion.
	aa=0.0
	av=0.0
	do i=1,ndat
		av=av+tinv(i) 
		aa=aa+ tinv(i)*tinv(i)
	end do
	av=av/ndat
	aa=aa/ndat
	a1=sqrt(aa-av*av)

	write(*,*)"residual_dev without errors=",a0
	write(*,*)"residual_dev with    errors=",a1
	end
!-----------------weighting the result-------------------
	subroutine weightout(bkw,resu,npar)
	dimension bkw(npar),resu(npar)
	ncol=npar
	do icol=1,ncol
		resu(icol)=resu(icol)*bkw(icol)
	end do
	end
!------------------------------------------------------------
	subroutine rsdchg(ndr,tim,nor,eqk,ndeq,sta,ndst,resu,ityp,ipar)
	include "surfpath.par"
	include "datatype.inc"
	dimension tim(3,ndr),nor(2,ndr),resu(npar)
	dimension nobk(nbk),segbk(nbk),azmbk(nbk),dmd(nbk,3)
	type (station)::sta(ndst)
	type (earthquake)::eqk(ndeq)
	type (partom)::ipar

	pp=3.1415926/180.0
	np=3
	kpar=ipar%ipar
	kdv=ipar%idv;
	if(kpar.ne.0)np=1
	fact=1.0
	do ia=1,ny
		il=(ia-1)*nx*np
		do io=1,nx
			ib=(ia-1)*nx+io
			ii=(ib-1)*np+1
			if (kdv.eq.2.0) then
				fact=-1.0/modva(ii,1)/modva(ii,1)
			end if
		    dmd(ib,1)=resu(ii)*fact

			if(np.lt.2)goto 100
			ii=(ib-1)*np+2
			dmd(ib,2)=resu(ii)
			ii=(ib-1)*np+3
			dmd(ib,3)=resu(ii)
100		end do
	end do

	do ir=1,ndr
		ie=nor(1,ir);elo=eqk(ie)%elo;ela=eqk(ie)%ela
		is=nor(2,ir);slo=sta(is)%slo;sla=sta(is)%sla
	    call bkofray(area,nx,ny,elo,ela,slo,sla,nbr,nobk,
     +		segbk,azmbk,re)
		dtinv=dtime0(kpar,dmd,nx,ny,nbr,nobk,segbk,azmbk)
		tim(3,ir)=tim(2,ir)+dtinv
	end do
	end
c---------------------------------------------------------------
	function dtime0(kpar,dmd,nx,ny,nbr,nobk,segbk,azmbk)
	dimension nobk(nbr),segbk(nbr),azmbk(nbr)
	real dmd(nx*ny,3)
	dt=0.0
	do ir=1,nbr
		nb=nobk(ir)
		ds=dmd(nb,1)
		da=dmd(nb,2)
		db=dmd(nb,3)
		seg=segbk(ir)
		az=azmbk(ir)
		if(kpar.eq.0)then
			dt=dt+seg*(ds+da*cos(2*az)+db*sin(2*az))
		elseif(kpar.eq.1)then
			dt=dt+seg*ds
		end if
	end do
	dtime0=dt
	end 

!------------------------------------------
	subroutine rglrz(ipar,wpn,wan,nrow,nz)
	include "surfpath.par"
	np=1
	if(ipar.eq.0)np=3   !if anisotropy included, parameter number is 3
	if(wpn.eq.0.0.and.wan.eq.0.0)then
		return
	elseif(wpn.eq.0.0)then
		ip1=2; ip2=np
		if(ip1.gt.ip2)return
	elseif(wan.eq.0.0)then
		ip1=1; ip2=1
	else
		ip1=1; ip2=np
	end if
	do 100 iy=1,ny
	do 100 ix=1,nx
	!how many points around the point(ix,iy) are in the net??
		near=0
		if((iy-1).ge.1)near=near+1
		if((ix-1).ge.1)near=near+1
		if((ix+1).le.nx)near=near+1
		if((iy+1).le.ny)near=near+1
	!add the smoothing constraint to the Jacobe Matrix??
		do ip=ip1,ip2
			nrow=nrow+1
			nnzr(nrow)=near+1
			wei=wan
			if(ip.eq.1)wei=wpn
			if(iy-1.ge.1)then
				ibk=(iy-1-1)*nx+ix
				call rgljcb(ibk,np,ip,wei,nz)
			end if
			if(ix-1.ge.1)then
				ibk=(iy-1)*nx+ix-1
				call rgljcb(ibk,np,ip,wei,nz)
			end if
			ibk=(iy-1)*nx+ix
			wxy=-wei*near
			call rgljcb(ibk,np,ip,wxy,nz)
			if(ix+1.le.nx)then
				ibk=(iy-1)*nx+ix+1
				call rgljcb(ibk,np,ip,wei,nz)
			end if
			if(iy+1.le.ny)then
				ibk=(iy+1-1)*nx+ix
				call rgljcb(ibk,np,ip,wei,nz)
			end if
		end do
100	end do
	end
