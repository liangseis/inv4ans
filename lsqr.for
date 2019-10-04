      SUBROUTINE LSQR(NOUT,itmax0,DAMP,M,N,B0,drmin,x,
     + nnzr,nzcl,nzvl,nry,nnz)
c      include "surfpath.par"      
      DIMENSION  U(M),B0(M),V(N),W(N),X(N),SE(N),resu(N)
	LOGICAL STP, INVLOC
       dimension nnzr(nry), nzcl(nnz)
       real nzvl(nnz)
        
c	NOUT---Output file logical unit
c	itmax--maximum iteration number
c	Damp==0.0
c	M------Number of rows
c	N---number of colcumn (number of parameters)
c	B0(M)---data vector (d=Gm)
c	drmin---exit critria (can be 0.0)
c	x----inversion result

      ATOL    = 0.0
      BTOL    = 0.0
      CONLIM  = 0.0
c      ITMAX   = 30

	DO IM=1,M
		U(IM)=B0(IM)
	END DO

      IF(NOUT.GT.0) WRITE(NOUT,1000) M,N,DAMP,ATOL,CONLIM,BTOL,ITMAX
c     initialization
      ZERO    = 0.0
      ONE     = 1.0
      CTOL    = ZERO
      IF(CONLIM.GT.ZERO) CTOL = ONE/CONLIM
      DAMPSQ  = DAMP*DAMP
      ANORM   = ZERO
      ACOND   = ZERO
      BBNORM  = ZERO
      DDNORM  = ZERO
      RES2    = ZERO
      XNORM   = ZERO
      XXNORM  = ZERO
      CS2     =-ONE
      SN2     = ZERO
      Z       = ZERO
      ITN     = 0
      ISTOP   = 0
      NSTOP   = 0
      DO 10 I = 1, N
      V(I)    = ZERO
      X(I)    = ZERO
      SE(I)   = ZERO
10    CONTINUE
      CALL NORMLZ(M,U,BETA)
      CALL APROD(2,M,N,V,U,nnzr,nzcl,nzvl,nry,nnz)
      CALL NORMLZ(N,V,ALFA)
      DO 20 I = 1,N
20    W(I)    = V(I)
      RHOBAR  = ALFA
      PHIBAR  = BETA
      BNORM   = BETA
      RNORM   = BETA
      ARNORM  = ALFA*BETA
      IF(ARNORM.LE.ZERO)       GO TO 800
      IF(NOUT.LE.0)            GO TO 100
      IF(DAMPSQ.LE.ZERO)       WRITE(NOUT, 1200)
      IF(DAMPSQ.GT.ZERO)       WRITE(NOUT, 1300)
      TEST1   = ONE
      TEST2   = ALFA/BETA
      WRITE(NOUT,1500) ITN,X(1),RNORM,TEST1,TEST2
      WRITE(NOUT,1600)
c     main iteration loop.
100   ITN     = ITN + 1
c     bidiagonalization     !
      AF      =-ALFA
      DO 30 I = 1,M
30    U(I)    = AF*U(I)
      CALL APROD(1,M,N,V,U,nnzr,nzcl,nzvl,nry,nnz)
      CALL NORMLZ(M,U,BETA)
      BBNORM  = BBNORM+ALFA*ALFA+BETA*BETA+DAMPSQ
      BT      =-BETA
      DO 40 I = 1,N
40    V(I)    = BT*V(I)
      CALL APROD(2,M,N,V,U,nnzr,nzcl,nzvl,nry,nnz)
      CALL NORMLZ(N,V,ALFA)
c     modified QR
      RHBAR2  = RHOBAR*RHOBAR+DAMPSQ
      RHBAR1  = SQRT(RHBAR2)
      CS1     = RHOBAR/RHBAR1
      SN1     = DAMP/RHBAR1
      PSI     = SN1*PHIBAR
      PHIBAR  = CS1*PHIBAR
      RHO     = SQRT(RHBAR2+BETA*BETA)
      CS      = RHBAR1/RHO
      SN      = BETA/RHO
      THETA   = SN*ALFA
      RHOBAR  =-CS*ALFA
      PHI     = CS*PHIBAR
      PHIBAR  = SN*PHIBAR
      TAU     = SN*PHI
c     update X, W and the standard error estimates.
      T1      = PHI/RHO
      T2      =-THETA/RHO
      T3      = ONE/RHO
      DO 50 I = 1, N
      T       = W(I)
      X(I)    = T1*T+X(I)
      W(I)    = T2*T+V(I)
      T3T     = T3*T
      T       = T3T*T3T
      SE(I)   = T+SE(I)
      DDNORM  = T+DDNORM
50    CONTINUE
      DELTA   = SN2*RHO
      GAMBAR  =-CS2*RHO
      RHS     = PHI-DELTA*Z
      ZBAR    = RHS/GAMBAR
      XNORM   = SQRT(XXNORM+ZBAR*ZBAR)
      GAMMA   = SQRT(GAMBAR*GAMBAR+THETA*THETA)
      CS2     = GAMBAR/GAMMA
      SN2     = THETA/GAMMA
      Z       = RHS/GAMMA
      XXNORM  = XXNORM + Z*Z
      ANORM   = SQRT(BBNORM)
      ACOND   = ANORM*SQRT(DDNORM)
      RES1    = PHIBAR*PHIBAR
      RES2    = RES2+PSI*PSI
      RNORM   = SQRT(RES1+RES2)
      ARNORM  = ALFA*ABS(TAU)
      TEST1   = RNORM/BNORM
      TEST2   = ARNORM/(ANORM*RNORM)
      TEST3   = ONE/ACOND
      T1      = TEST1/(ONE+ANORM*XNORM/BNORM)
      RTOL    = BTOL+ATOL*ANORM*XNORM/BNORM
      T3      = ONE+TEST3
      T2      = ONE+TEST2
      T1      = ONE+T1

c	CALL limit(ITN,ludv,updv,N,invloc,X)
c	CALL stopit(ITN,NOUT,M,N,drmin,X,b0,stp,resu)
      WRITE(NOUT,1250) ITN,ACOND,T1,T2,T3,RNORM
c      write(6,1250) ITN,ACOND,T1,T2,T3,RNORM
1250  FORMAT(I4,4F12.6,F10.3,f10.4)

	IF(STP)write(*,*)"ok to stop"
      IF(ITN.LT.10)             GO TO 100
      IF(ITN.GE.ITMAX)          ISTOP = 7
      IF(T3.LE.ONE)             ISTOP = 6
      IF(T2.LE.ONE)             ISTOP = 5
      IF(T1.LE.ONE)             ISTOP = 4
      IF(TEST3.LE.CTOL)         ISTOP = 3
      IF(TEST2.LE.ATOL)         ISTOP = 2
      IF(TEST1.LE.RTOL)         ISTOP = 1
      IF(NOUT.LE.0)             GO TO 600
      IF(M.LE.40.OR.N.LE.40)    GO TO 400
      IF(ITN.LE.10)             GO TO 400
      IF(ITN.GE.ITMAX-10)       GO TO 400
      IF(MOD(ITN,10).EQ.0)      GO TO 400
      IF(TEST3.LE.2.0*CTOL)     GO TO 400
      IF(TEST2.LE.10.0*ATOL)    GO TO 400
      IF(TEST1.LE.10.0*RTOL)    GO TO 400
      GO TO 600
400   CONTINUE
c     WRITE(NOUT,1500) ITN,X(1),RNORM,TEST1,TEST2,ANORM,ACOND
      IF(MOD(ITN,10).EQ.0)      WRITE(NOUT, 1600)
600   IF(ISTOP.EQ.0)            NSTOP = 0
      IF(ISTOP.EQ.0)            GO TO 100
      NCONV   = 1
      NSTOP   = NSTOP+1
      IF(NSTOP.LT.NCONV.AND.ITN.LT.ITMAX)  ISTOP = 0
      IF(ISTOP.EQ.0)            GO TO 100
c     end of iteration loop
      T       = ONE
      IF(M.GT.N)                T = M - N
      IF(DAMPSQ.GT.ZERO)        T = M
      T       = RNORM/SQRT(T)
      DO 70 I = 1, N
      SE(I)   = T*SQRT(SE(I))
70    CONTINUE
800   IF(NOUT.LE.0)             GO TO 900
      WRITE(NOUT,1900)          ITN,ISTOP
      IF(ISTOP.EQ.0)            WRITE(NOUT,2000)
      IF(ISTOP.EQ.1)            WRITE(NOUT,2100)
      IF(ISTOP.EQ.2)            WRITE(NOUT,2200)
      IF(ISTOP.EQ.3)            WRITE(NOUT,2300)
      IF(ISTOP.EQ.4)            WRITE(NOUT,2400)
      IF(ISTOP.EQ.5)            WRITE(NOUT,2500)
      IF(ISTOP.EQ.6)            WRITE(NOUT,2600)
      IF(ISTOP.EQ.7)            WRITE(NOUT,2700)
 1000 FORMAT(// 25X,"LSQR   --   Least-squares solution of  A*X = B"
     &   // 25X,"The matrix  A  has ", I6,"rows and ", I6," cols"
     &    / 25X,"The damping parameter is DAMP = ", 1PE10.2
     &   // 25X,"ATOL  = ", 1PE10.2, 10X, " CONLIM = ", 1PE10.2
     &    / 25X,"BTOL  = ", 1PE10.2, 10X, " ITMAX  = ", I10)
 1200 FORMAT(// 3X, 3HITN, 9X, 4HX(1), 14X, 8HFUNCTION, 7X,
     &   "COMPATIBLE INCOMPATIBLE    NORM(A)  COND(A)"/)
 1300 FORMAT(// 3X, 3HITN, 9X, 4HX(1), 14X, 8HFUNCTION, 7X,
     &   "COMPATIBLE INCOMPATIBLE NORM(ABAR) COND(ABAR)"/)
 1500 FORMAT(I6, 1PE20.10, 1PE19.10, 1P2E13.3, 1P2E11.2)
 1600 FORMAT(1X)
 1900 FORMAT(/" No. of iterations =",I6,8X," Stopping condition =",I3)
 2000 FORMAT(/" The exact solution is  X = 0.")
 2100 FORMAT(/" A*X - B  is small enough, given  ATOL, BTOL")
 2200 FORMAT(/" The least-sqrs solu is good enough, given  ATOL")
 2300 FORMAT(/" The estimate of  COND(ABAR)  has exceeded  CONLIM")
 2400 FORMAT(/" A*X - B  is small enough for this machine")
 2500 FORMAT(/" The least-sqrs solu is good enough for this machine")
 2600 FORMAT(/" COND(ABAR)  seems to be too large for this machine")
 2700 FORMAT(/" The iteration limit has been reached")
 900  RETURN
      END
!------------------------
      SUBROUTINE APROD(MODE,M,N,X,Y,nnzr,nzcl,nzvl,nry,nnz)
c      INCLUDE "PARAM" 
c      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
c     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
       dimension nnzr(nry), nzcl(nnz)
       real nzvl(nnz)

c      include "surfpath.par"

      DIMENSION  X(N),Y(M)
      ZERO   = 0.0
      L2     = 0
      IF(MODE.NE.1)  GO TO 4
c     MODE = 1 -- SET  Y = Y+A*X.
      DO 1 I = 1,M
      SUM    = ZERO
      L1     = L2+1
      L2     = L2+nnzr(I)
      DO 2 L = L1,L2
      J      = nzcl(L)
      SUM    = SUM+nzvl(L)*X(J)
2     CONTINUE
      Y(I)   = Y(I)+SUM
1     CONTINUE
      RETURN
4     CONTINUE
c     MODE = 2 -- SET  X = X+AT*Y
      DO 5 I = 1,M
      YI     = Y(I)
      L1     = L2+1
      L2     = L2+nnzr(I)
      DO 6 L = L1,L2
      J      = nzcl(L)
      X(J)   = X(J)+nzvl(L)*YI
6     CONTINUE
5     CONTINUE
      RETURN
      END
!---------------
      SUBROUTINE NORMLZ(N,X,S)
      DIMENSION  X(N)
c     normalizes vector X
      EPS    = 1.0E-10
      S      = 0.0
      DO 1 I = 1,N
      XI     = X(I)
1     S      = S+XI*XI
      S      = SQRT(S)
      SS     = 0.0
      IF(S.GT.EPS) SS = 1.0/S
      DO 2 I = 1,N
2     X(I)   = X(I)*SS
      RETURN
      END
c------------------------------------------------------------------
c$$$	subroutine stopit(its,lulg,ndata,ngl,drmin,solu,b0,stp,resu)
c$$$	include "surfpath.par"
c$$$	logical stp
c$$$	dimension resu(NGL)
c$$$
c$$$	stp=.false.
c$$$      dax=0.0
c$$$	dr=0.0
c$$$	NELEA=0
c$$$	DO  I=1,NDATA
c$$$		L=nnzr(I)
c$$$		xx=0.0
c$$$		DO  L1=1,L
c$$$			NELEA=NELEA+1
c$$$			L2=nzcl(NELEA)
c$$$			xx=xx+nzvl(NELEA)*solu(L2)
c$$$		end do
c$$$		rr=b0(i)-xx
c$$$		dax=dax+xx*xx
c$$$		dr=dr+rr*rr
c$$$      end do
c$$$	dax=sqrt(dax/(ndata-2));  dr=sqrt(dr/(ndata-2))
c$$$
c$$$	de=dax/dr
c$$$	if(dr.gt.drmin)then
c$$$	    stp=.true.
c$$$	else
c$$$           do ik=1,NGL
c$$$              resu(ik)=solu(ik)
c$$$           end do
c$$$           drmin=dr
c$$$	end if
c$$$
c$$$c	write(lulg,'(i4,3f10.6)')its,dax,dr,de
c$$$	end
