! -*- F90 -*-
      subroutine hknevolve(x,Q,f) 
      implicit none 
      integer nq,nx,nd,nff,nset, nhess
      PARAMETER (NQ=33, NX=117, ND=7, NFF=7, nhess=0)
      include 'parmsetup.inc'
      character*16 name(nmxset)
      character*512 setpath
      integer nmem(nmxset),ndef(nmxset),mmem,mem
      common/NAME/name,nmem,ndef,mmem
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
      CHARACTER*80 LINE
      double precision pdf(-6:6),x,q,q2
      double precision QG(NQ),XG(NX),PDFG(NX,NQ,ND,0:nhess),DNPDF(-4:4) &
     & ,BXG(NX,NQ,ND,0:nhess), CXG(NX,NQ,ND,0:nhess), DXG(NX,NQ,ND,0:nhess)
      double precision PDFJ1(ND), PDFJ2(ND)
      double precision T,DX
      integer iset,imem 
      integer i,j,k,n
      integer iserch
      real*8 f(-6:6) 
      real*8 alfas 
      real*8 Eorder,Q2fit 
      DATA QG / &
     &  1.000000D+00, 1.467799D+00, 2.154435D+00, &
     &  3.162278D+00, 4.641589D+00, 6.812921D+00, &
     &  1.000000D+01, 1.467799D+01, 2.154435D+01, &
     &  3.162278D+01, 4.641589D+01, 6.812921D+01, &
     &  1.000000D+02, 1.778279D+02, 3.162278D+02, 5.623413D+02, &
     &  1.000000D+03, 1.778279D+03, 3.162278D+03, 5.623413D+03, &
     &  1.000000D+04, 1.778279D+04, 3.162278D+04, 5.623413D+04, &
     &  1.000000D+05, 1.778279D+05, 3.162278D+05, 5.623413D+05, &
     &  1.000000D+06, 4.641589D+06, & 
     &  1.000000D+07, 4.641589D+07, & 
     &  1.000000D+08  /

      DATA XG / &
     &  1.000000D-09, 1.333521D-09, 1.778279D-09, 2.371374D-09, &
     &  3.162278D-09, 4.216965D-09, 5.623413D-09, 7.498942D-09, &
     &  1.000000D-08, 1.333521D-08, 1.778279D-08, 2.371374D-08, &
     &  3.162278D-08, 4.216965D-08, 5.623413D-08, 7.498942D-08, &
     &  1.000000D-07, 1.333521D-07, 1.778279D-07, 2.371374D-07, &
     &  3.162278D-07, 4.216965D-07, 5.623413D-07, 7.498942D-07, &
     &  1.000000D-06, 1.333521D-06, 1.778279D-06, 2.371374D-06, &
     &  3.162278D-06, 4.216965D-06, 5.623413D-06, 7.498942D-06, &
     &  1.000000D-05, 1.333521D-05, 1.778279D-05, 2.371374D-05, &
     &  3.162278D-05, 4.216965D-05, 5.623413D-05, 7.498942D-05, &
     &  1.000000D-04, 1.333521D-04, 1.778279D-04, 2.371374D-04, &
     &  3.162278D-04, 4.216965D-04, 5.623413D-04, 7.498942D-04, &
     &  1.000000D-03, 1.154782D-03, 1.333521D-03, 1.539927D-03, &
     &  1.778279D-03, 2.053525D-03, 2.371374D-03, 2.738420D-03, &
     &  3.162278D-03, 3.651741D-03, 4.216965D-03, 4.869675D-03, &
     &  5.623413D-03, 6.493816D-03, 7.498942D-03, 8.659643D-03, &
     &  1.000000D-02, 1.154782D-02, 1.333521D-02, 1.539927D-02, &
     &  1.778279D-02, 2.053525D-02, 2.371374D-02, 2.738420D-02, &
     &  3.162278D-02, 3.651741D-02, 4.216965D-02, 4.869675D-02, &
     &  5.623413D-02, 6.493816D-02, 7.498942D-02, 8.659643D-02, &
     &  1.000000D-1, 1.250000D-1, 1.500000D-1, 1.750000D-1, &
     &  2.000000D-1, 2.250000D-1, 2.500000D-1, 2.750000D-1, &
     &  3.000000D-1, 3.250000D-1, 3.500000D-1, 3.750000D-1, &
     &  4.000000D-1, 4.250000D-1, 4.500000D-1, 4.750000D-1, &
     &  5.000000D-1, 5.250000D-1, 5.500000D-1, 5.750000D-1, &
     &  6.000000D-1, 6.250000D-1, 6.500000D-1, 6.750000D-1, &
     &  7.000000D-1, 7.250000D-1, 7.500000D-1, 7.750000D-1, &
     &  8.000000D-1, 8.250000D-1, 8.500000D-1, 8.750000D-1, &
     &  9.000000D-1, 9.250000D-1, 9.500000D-1, 9.750000D-1, &
     &  1.000000D+0 /
     
      save 
      
      call getnset(iset)
      call getnmem(iset,imem)
     
      DO I=-4,4
        DNPDF(I)=0.D0
      END DO
  
      q2 = q*q
! CHECK X AND Q2 VALUES.
      IF((X.LT.1.D-9).OR.(X.GT.1.D0)) THEN
        WRITE(*,1030) X
 1030   FORMAT (' ','FF WARNING: OUT OF RANGE --> X =', 1PE12.3)
        STOP
      ENDIF
      IF((Q2.LT.1.D0).OR.(Q2.GT.1.D8)) THEN
        WRITE(*,1040) Q2
 1040   FORMAT (' ','FF WARNING: OUT OF RANGE --> Q2 =', 1PE12.3)
        STOP
      ENDIF

! INTERPOLATION.
! X: CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION.
      J=ISERCH(NQ,QG,Q2)
      IF(J.EQ.NQ) J=NQ-1
      K=ISERCH(NX,XG,X)
      DO I=1,ND
        DX=X-XG(K)
        PDFJ1(I)=PDFG(K,J,I,0) &
     &       +DX*(BXG(K,J,I,0)+DX*(CXG(K,J,I,0)+DX*DXG(K,J,I,0)))
        PDFJ2(I)=PDFG(K,J+1,I,0) &
     &       +DX*(BXG(K,J+1,I,0)+DX*(CXG(K,J+1,I,0)+DX*DXG(K,J+1,I,0)))

      ENDDO

! -- Nuclear PDF functions --
      T=(DLOG(Q2)-DLOG(QG(J)))/(DLOG(QG(J+1))-DLOG(QG(J)))
      f(0)=(1.D0-T)*PDFJ1(1)+T*PDFJ2(1)     ! g
      f(1)=(1.D0-T)*PDFJ1(3)+T*PDFJ2(3)     ! d
      f(2)=(1.D0-T)*PDFJ1(2)+T*PDFJ2(2)     ! u
      f(-1)=(1.D0-T)*PDFJ1(5)+T*PDFJ2(5)     ! db
      f(-2)=(1.D0-T)*PDFJ1(4)+T*PDFJ2(4)     ! ub
      f(-3)=(1.D0-T)*PDFJ1(6)+T*PDFJ2(6)     ! sb
      f(4)=(1.D0-T)*PDFJ1(7)+T*PDFJ2(7) !  c
      f(3)=f(-3) ! s=sb
      f(-4)=f(4) ! cb=c
      return 

!                                                                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      
      entry hkngetgrid(nset,ngridx,ngridq,gridx,gridq)
     
      ngridx=NX
      do jx=1,ngridx
          gridx(jx)=XG(jx)
      enddo
      ngridq=NQ
      do jq=1,ngridq
          gridq(jq)=QG(jq)
      enddo
       
      return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      entry hknread(nset) 
  !dummy read to get to the End: (stream 1 is still open)

      read(1,*)nmem(nset),ndef(nset)

      do n=0,nmem(nset)
        DO J=1,NQ
          DO K=1,NX-1
            READ(1,'(a)') line
          ENDDO
        ENDDO
      enddo
 
      return 
!                                                                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      entry hknalfa(alfas,Q)
!        call alphamrs(4,alfas,q)
        call alphahkn(q,alfas)
      return 
!                                                                          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      entry hkninit(nset,Eorder,Q2fit) 
      return 
!                                                                       
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      entry hknpdf(mem) 

    call getnset(iset)
    call setnmem(iset,mem)

  ! here is the real read!!

  ! have to reopen stream 1 
  call getsetpath(setpath)
  open(1,file=setpath(1:len_trim(setpath)),action='READ')

  line = ''
  do while (line(2:11).ne.'Evolution:') 
     read(1,'(a)'),line
  enddo
  read(1,'(a)'),line
  read(1,'(a)'),line

  read(1,*)nmem(iset),ndef(iset)
  ! - dummy read up to the member requested
  do i=0,mem-1
     do j=1,nq
        do k=1,nx-1
           read(1,'(a)')line
        enddo
     enddo
  enddo
  !   Now read in the grids from the grid file.
      DO J=1,NQ
        DO K=1,NX-1
          READ(1,1025) (PDFG(K,J,I,0), I=1,NFF)
!         print *,
        ENDDO
      ENDDO
      DO I=1,ND
        DO J=1,NQ
          PDFG(NX,J,I,0)=0.D0 ! x=1 NPDF=0.D0
          CALL LSPLINE(NX,XG,PDFG,BXG,CXG,DXG,ISET,I,J,0)
        ENDDO
      ENDDO
      
      close(1) 

 1025 FORMAT(1X,7(1PE14.6))

      return 
!                                                                       
      END                                           
! ---------------------------------------------------------------------
      SUBROUTINE LSPLINE(N,X,Y,B,C,D,ISET,I,J,nmem)
! ---------------------------------------------------------------------
! CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
! INTERPOLATION SUBROUTINES ARE TAKEN FROM
! G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
! COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NQ=33, NX=117, ND=7, nhess=12)
      DIMENSION Y(NX,NQ,ND,0:nhess),B(NX,NQ,ND,0:nhess),C(NX,NQ,ND,0:nhess),D(NX,NQ,ND,0:nhess) &
     &         ,X(NX) 
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1,J,I,nmem)=X(2)-X(1)
      C(2,J,I,nmem)=(Y(2,J,I,nmem)-Y(1,J,I,nmem))/D(1,J,I,nmem)
      DO 210 K=2,NM1
        D(K,J,I,nmem)=X(K+1)-X(K)
        B(K,J,I,nmem)=2.0D0*(D(K-1,J,I,nmem)+D(K,J,I,nmem))
        C(K+1,J,I,nmem)=(Y(K+1,J,I,nmem)-Y(K,J,I,nmem))/D(K,J,I,nmem)
        C(K,J,I,nmem)=C(K+1,J,I,nmem)-C(K,J,I,nmem)
  210 CONTINUE
      B(1,J,I,nmem)=-D(1,J,I,nmem)
      B(N,J,I,nmem)=-D(N-1,J,I,nmem)
      C(1,J,I,nmem)=0.0D0
      C(N,J,I,nmem)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1,J,I,nmem)=C(3,J,I,nmem)/(X(4)-X(2))-C(2,J,I,nmem)/(X(3)-X(1))
      C(N,J,I,nmem)=C(N-1,J,I,nmem)/(X(N)-X(N-2))-C(N-2,J,I,nmem)/(X(N-1)-X(N-3))
      C(1,J,I,nmem)=C(1,J,I,nmem)*D(1,J,I,nmem)**2.0D0/(X(4)-X(1))
      C(N,J,I,nmem)=-C(N,J,I,nmem)*D(N-1,J,I,nmem)**2.0D0/(X(N)-X(N-3))
  215 CONTINUE
      DO 220 K=2,N
        T=D(K-1,J,I,nmem)/B(K-1,J,I,nmem)
        B(K,J,I,nmem)=B(K,J,I,nmem)-T*D(K-1,J,I,nmem)
        C(K,J,I,nmem)=C(K,J,I,nmem)-T*C(K-1,J,I,nmem)
  220 CONTINUE
      C(N,J,I,nmem)=C(N,J,I,nmem)/B(N,J,I,nmem)
      DO 230 IB=1,NM1
        K=N-IB
        C(K,J,I,nmem)=(C(K,J,I,nmem)-D(K,J,I,nmem)*C(K+1,J,I,nmem))/B(K,J,I,nmem)
  230 CONTINUE
      B(N,J,I,nmem)=(Y(N,J,I,nmem)-Y(NM1,J,I,nmem))/D(NM1,J,I,nmem) &
     &        +D(NM1,J,I,nmem)*(C(NM1,J,I,nmem)+2.0D0*C(N,J,I,nmem))
      DO 240 K=1,NM1 
        B(K,J,I,nmem)=(Y(K+1,J,I,nmem)-Y(K,J,I,nmem))/D(K,J,I,nmem) &
     &          -D(K,J,I,nmem)*(C(K+1,J,I,nmem)+2.0D0*C(K,J,I,nmem))
        D(K,J,I,nmem)=(C(K+1,J,I,nmem)-C(K,J,I,nmem))/D(K,J,I,nmem)
        C(K,J,I,nmem)=3.0D0*C(K,J,I,nmem)
  240 CONTINUE
      C(N,J,I,nmem)=3.0D0*C(N,J,I,nmem)
      D(N,J,I,nmem)=D(N-1,J,I,nmem)
      RETURN
  250 CONTINUE
      B(1,J,I,nmem)=(Y(2,J,I,nmem)-Y(1,J,I,nmem))/(X(2)-X(1))
      C(1,J,I,nmem)=0.0D0
      D(1,J,I,nmem)=0.0D0
      B(2,J,I,nmem)=B(1,J,I,nmem)
      C(2,J,I,nmem)=0.0D0
      D(2,J,I,nmem)=0.0D0
      RETURN
      END
! ---------------------------------------------------------------------
      INTEGER FUNCTION ISERCH(N,X,Y)
! ---------------------------------------------------------------------
! THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
! X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(117)

      MIN=1
      MAX=N+1

   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10

      ISERCH=MIN

      RETURN
      END
! *********************************************************************
! THE END OF THE PROGRAM.
! *********************************************************************
! ---------------------------------------------------------------------
!  IN: Q2=Q^2 [GeV^2], IORDER=1:LO, 2:NLO 
! OUT: Alpha_s
!
      SUBROUTINE alphahkn(Q,ALPHA_S)
! ---------------------------------------------------------------------
! RUNNING COUPLING CONSTANTS.
      IMPLICIT REAL*8(A-H,O-Z)
      DATA DLAML,DLAMN/0.174D0, 0.3D0/
!      DATA THRE4,THRE5/1.35D0, 4.3D0/
      DATA THRE4,THRE5/1.35D0, 1.D+6/

      call getnset(nset) 
      call GetOrderAsM(nset,iord)
      q2 = q*q 
      iorder = iord + 1    
 
      PI=4.D0*DATAN(1.D0)
      CTHRE=THRE4*THRE4
      BTHRE=THRE5*THRE5
      CF=4.D0/3.D0
      CG=3.D0
      TR=1.D0/2.D0

! Changing the number of the quark flavor at heavy quark mass threshold
      Q2thr=Q2 
      IF(Q2thr.LT.CTHRE) F=3.D0  
      IF((Q2thr.GE.CTHRE).AND.(Q2thr.LT.BTHRE)) F=4.D0
      IF(Q2thr.GE.BTHRE) F=5.D0                      

      B0=11.D0/3.D0*CG-4.D0/3.D0*TR*F
      B1=34.D0/3.D0*CG*CG-10.D0/3.D0*CG*F-2.D0*CF*F

! Changing Lambda_QCD for connecting alpah_s at the threshold 
      IF(Q2thr.LT.CTHRE) then      ! Lambda_QCD (4to3)
        DLAMF=DLAML*(DSQRT(CTHRE)/DLAML)**(2.D0/3.D0/B0)
        DLAMFN=DLAMN*(DSQRT(CTHRE)/DLAMN)**(2.D0/27.D0) &
     &        *DLOG(CTHRE/(DLAMN*DLAMN))**(107.D0/2025.D0)

      ELSE IF((Q2thr.GE.CTHRE).AND.(Q2thr.LT.BTHRE)) then 
        DLAMF=DLAML                 
        DLAMFN=DLAMN

      ELSE IF(Q2thr.GE.BTHRE) Then ! Lambda_QCD (4to5) 
        DLAMF=DLAML*(DSQRT(BTHRE)/DLAML)**(-2.D0/3.D0/B0)
        DLAMFN=DLAMN*(DLAMN/DSQRT(BTHRE))**(2.D0/23.D0) &
     &        *DLOG(BTHRE/(DLAMN*DLAMN))**(-963.D0/13225.D0)
      END IF

! Calculating alpha_s(Q^2)
      IF(IORDER.EQ.2) DLAMF=DLAMFN
      DLNLAM=DLOG(DLAMF*DLAMF)
      DLNQ2=DLOG(Q2)-DLNLAM

      ALPHA=4.D0*PI/B0/DLNQ2    ! LO
      IF(IORDER.EQ.2) THEN      ! NLO
        ALPHA=ALPHA*(1.D0-B1*DLOG(DLNQ2)/(B0*B0*DLNQ2))
      END IF

      ALPHA_S=ALPHA

      RETURN
      END
! ---------------------------------------------------------------------
