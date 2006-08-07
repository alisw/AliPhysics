      subroutine ZEUSevolve(x,Q,pdf)
      implicit double precision (a-h,o-z)
      include 'parmsetup.inc'
      character*64 gridname
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem
      common/NAME/name,nmem,ndef,mem
      integer nset,iset,isetlast
      data isetlast/-1/
      integer Eorder
      real*8 mc,mc2,mb,mb2,mt,mt2
      real*8 f(-6:6),pdf(-6:6)
      integer qnerr
      parameter(nstartp=7)
      DIMENSION QSP(NSTARTP)
      DATA QSP/10.,20.,30.,40.,50.,80.,100./
      real*8 xval(45)
      logical HEAVY,VFN
      real*8 pwgt(20)
      save
*
      call getnset(iset)
c      print *,'iset=',iset,' now calling get_pdfqcd'
      if(iset.ne.isetlast) then
        call get_pdfqcd(iset)
        isetlast = iset
      endif
*
      Q2=Q*Q
      UCENT=QPDFXQ('UPVAL',X,Q2,IFAIL)
      DCENT=QPDFXQ('DNVAL',X,Q2,IFAIL)
      GCENT=QPDFXQ('GLUON',X,Q2,IFAIL)
      UBCEN=QPDFXQ('UB',X,Q2,IFAIL)
      DBCEN=QPDFXQ('DB',X,Q2,IFAIL)
      STCEN=QPDFXQ('SB',X,Q2,IFAIL)
      IF(Q2.GE.Q2C)THEN
        CHCEN=QPDFXQ('CB',X,Q2,IFAIL)
      ELSE
        CHCEN=0.0
      ENDIF
      IF(Q2.GE.Q2B)THEN
        BTCEN=QPDFXQ('BB',X,Q2,IFAIL)
      ELSE
        BTCEN=0.0
      ENDIF
      pdf(1)=dcent+dbcen
      pdf(2)=ucent+ubcen
      pdf(3)=stcen
      pdf(4)=chcen
      pdf(5)=btcen
      pdf(6)=0d0
      pdf(0)=gcent
      pdf(-1)=dbcen
      pdf(-2)=ubcen
      pdf(-3)=stcen
      pdf(-4)=chcen
      pdf(-5)=bbtcen
      pdf(-6)=0d0
*
      return
*
c===========================================================================
      entry ZEUSalfa(alfas,Q)
      Q2=Q*Q
      nf=6
      if (Q2.lt.mt2) nf=5
      if (Q2.lt.mb2) nf=4
      if (Q2.lt.mc2) nf=3
      alfas=QALFAS(Q2,Qlam,nf,iflag)
c      print *,q2,alfas,qlam,nf,iflag
      return
*
c===========================================================================
      entry ZEUSread(nset)
         read(1,*) gridname,nx,xmin,xmax,nq,qmin,qmax
      return
*
c===========================================================================
      entry ZEUSinit(nset,Eorder,Q2fit)
c      print *,name(nset)
      if(name(nset).eq.'QCDNUM_ZEUS_TR') then
         HEAVY = .FALSE.
	 VFN = .TRUE.
      else if(name(nset).eq.'QCDNUM_ZEUS_FF') then
         HEAVY = .TRUE.
	 VFN = .FALSE.
      else if(name(nset).eq.'QCDNUM_ZEUS_ZM') then
         HEAVY = .FALSE.
	 VFN = .FALSE.
      else
        print *,'name/scheme not recognized'
	stop 1
      endif
c--try 3 way logic ffn/zm-vfn/rt-vfn
      IF(HEAVY)THEN
      IVFN=1
      ELSE
      IVFN=0
      ENDIF
      IF(VFN)THEN
      IVFN=IVFN+2
      ELSE
      IVFN=IVFN
      ENDIF
C IVFN=0 IS ZM-VFN, 1 IS FFN,2 IS RT-VFN, 3 IS NOT ALLOWED  

      IF(IVFN.EQ.3)THEN
      WRITE(*,*)'IVFN=3 SO STOP',IVFN
      STOP
      ENDIF    

      
      
c--qcdnum initialisation
      CALL QNINIT
c--se thresholds
      Q0=Q2fit
      ZM=91.187D0
      ZM2=ZM*ZM
      ALPHAS=QNALFA(ZM2)

      call getQmassM(nset,4,mc)
      mc2=mc*mc
      call getQmassM(nset,5,mb)
      mb2=mb*mb
      call getQmassM(nset,6,mt)
      mt2=mt*mt
    
c      Q2C=1.8225
      Q2C=mc2
c      Q2B=18.49
      Q2B=mb2 
c      print  *,q2c,q2b    

      IF (Q0.LT.Q2C) THEN
        NACT=3
      ELSE
        NACT=4
      ENDIF
c--this merely defines nact where we startevolution
c--namely at q0
      IF (HEAVY) NACT=3

      CALL QNRSET('MCSTF',SQRT(Q2C))
      CALL QNRSET('MBSTF',SQRT(Q2B))
      CALL QNRSET('MCALF',SQRT(Q2C))
      CALL QNRSET('MBALF',SQRT(Q2B))
    
      IF (HEAVY) THEN
        CALL QTHRES(1D10,2D10)
c        CALL QTHRES(1D6,2D6)
      ELSE
        CALL QTHRES(Q2C,Q2B)
      ENDIF
 
      DO I=1,NSTARTP
        CALL GRQINP(QSP(I),1)
      ENDDO
      CALL GRQINP(Q0,1)
      CALL GRQINP(Q2C,1)
      CALL GRQINP(Q2B,1)
c      qcdnum grid not my grid

c      CALL GRXLIM(120,97D-8)
      CALL GRXLIM(nx,xmin,xmax)

c      CALL GRQLIM(61,29D-2,200D3) 
      CALL GRQLIM(nq,qmin,qmax) 

C--   Get final grid definitions and grid indices of Q0, Q2C and Q2B

      CALL GRGIVE(NXGRI,XMI,XMA,NQGRI,QMI,QMA)
c      WRITE(*,*)'NX,XL,XH,NQ,QL,QH',NXGRI,XMI,XMA,NQGRI,QMI,QMA
      IQ0 = IQFROMQ(Q0)
      IQC = IQFROMQ(Q2C)
      IQB = IQFROMQ(Q2B)
C--   Allow for heavy weights

      IF (HEAVY) THEN
        CALL QNLSET('WTF2C',.TRUE.)
        CALL QNLSET('WTF2B',.TRUE.)
        CALL QNLSET('CLOWQ',.FALSE.)
        CALL QNLSET('WTFLC',.TRUE.)
        CALL QNLSET('WTFLB',.TRUE.)
      ENDIF

C--   Compute weights and dump, or read in
c
c      IF (READIN) THEN 
c        OPEN(UNIT=24,FILE='weights.dat',FORM='UNFORMATTED',
c     .                                  STATUS='UNKNOWN')
c        CALL QNREAD(24,ISTOP,IERR)
c      ELSE
c        CALL QNFILW(0,0)
c        IF (HEAVY) THEN
c          OPEN(UNIT=24,FILE='weights.dat',FORM='UNFORMATTED',
c     .                                    STATUS='UNKNOWN')
c          CALL QNDUMP(24)
c        ENDIF
c      ENDIF


      if (index(gridname,'none').eq.1) then
         call qnfilw(0,0)
      else
         qnerr=-1
         open(unit=2,status='old',file=gridname,
     .        form='unformatted',err=1)
         call QNREAD(2,1,qnerr)
 1       close(2)
         if (qnerr.ne.0) then
            write(*,*) 'Grid file problem: ',gridname
            if (qnerr.lt.0) then 
               write(*,*) 'Grid file does not exist'
               write(*,*) 'Calculating and creating grid file'
               call qnfilw(0,0)
               open(unit=2,status='unknown',file=gridname,
     .              form='unformatted')
               call QNDUMP(2)
               close(2)
            else
               write(*,*) 'Existing grid file is inconsistent'
               if (qnerr.eq.1)
     .              write(*,*) 'Defined grid different'
               if (qnerr.eq.2)
     .              write(*,*) 'Heavy quark weight table different'
               if (qnerr.eq.3)
     .              write(*,*) 'Charm mass different'
               if (qnerr.eq.4)
     .              write(*,*) 'Bottom mass different'
               stop
            endif
         endif
      endif

C--   Apply cuts to grid
c--taking away the s cut at 600d0
      CALL GRCUTS(-1D0,-1D0,-1D0,-1D0)




C--   Choose renormalisation and factorisation scales

      CALL QNRSET('AAAR2',1D0)  ! renormalisation
      CALL QNRSET('BBBR2',0D0)
      CALL QNRSET('AAM2L',1D0)  ! factorisation (light)
      CALL QNRSET('BBM2L',0D0)
      CALL QNRSET('AAM2H',1D0)  ! factorisation (heavy)
      CALL QNRSET('BBM2H',0D0)
 
c       ZM=91.187D0
      imem=0
c      print *,imem
c -- only need call to listPDF here to get alphas
      call listPDF(nset,imem,xval)
c      print *,xval
        AS=XVAL(1)
c        AS=0.118d0
      CALL QNRSET('ALFQ0',ZM*ZM)
      CALL QNRSET('ALFAS',AS)

c      ZM2=ZM*ZM
      ALPHAS=QNALFA(ZM2)
c      WRITE(*,*)'ALPHAS AT Mz2',ALPHAS

C--   Book non-singlet distributions

      CALL QNBOOK(2,'UPLUS')
      CALL QNBOOK(3,'DPLUS')
      CALL QNBOOK(4,'SPLUS')
      CALL QNBOOK(5,'CPLUS')
      CALL QNBOOK(6,'BPLUS')
      CALL QNBOOK(7,'UPVAL')
      CALL QNBOOK(8,'DNVAL')
 
C--   Book linear combinations for proton for f = 3,4,5 flavours

c--define some quark pdfs
         CALL dVZERO(PWGT,20)        
        PWGT(2) = 0.5

        PWGT(7) = -0.5

        PWGT(1) = 0.5/3.
        CALL QNLINC(17,'UB',3,PWGT)
        PWGT(1) = 0.5/4.
        CALL QNLINC(17,'UB',4,PWGT)
        PWGT(1) = 0.5/5.
        CALL QNLINC(17,'UB',5,PWGT) 
        CALL dVZERO(PWGT,20) 

        PWGT(4) = 0.5
        PWGT(1) = 0.5/3.
        CALL QNLINC(18,'SB',3,PWGT)
        PWGT(1) = 0.5/4.
        CALL QNLINC(18,'SB',4,PWGT)
        PWGT(1) = 0.5/5.
        CALL QNLINC(18,'SB',5,PWGT)
         CALL dVZERO(PWGT,20)        
        CALL QNLINC(19,'CB',3,PWGT)
        PWGT(5) = 0.5
        PWGT(1) = 0.5/4.

        CALL QNLINC(19,'CB',4,PWGT)
        PWGT(1) = 0.5/5.
        CALL QNLINC(19,'CB',5,PWGT) 
        CALL dVZERO(PWGT,20) 
        PWGT(3) = 0.5

        PWGT(8) = -0.5

        PWGT(1) = 0.5/3.
        CALL QNLINC(20,'DB',3,PWGT)
        PWGT(1) = 0.5/4.
        CALL QNLINC(20,'DB',4,PWGT)
        PWGT(1) = 0.5/5.
        CALL QNLINC(20,'DB',5,PWGT)
         CALL dVZERO(PWGT,20)        
        CALL QNLINC(21,'BB',3,PWGT)
       CALL QNLINC(21,'BB',4,PWGT)  
        PWGT(6) = 0.5
  
        PWGT(1) = 0.5/5.
        CALL QNLINC(21,'BB',5,PWGT) 
c---

      return
*
c==========================================================================
      entry ZEUSpdf(nset)
c      ZM = 91.187d0
c      zm2 = zm*zm

c      ALPHAS=QNALFA(ZM2)
      
      
c      imem=mem
      call getnmem(nset,imem)    
      call listPDF(nset,imem,xval)
c      print *,nset,imem
c      print *,xval
c      print *,imem,xval
      
      UA=XVAL(3)
      UB=XVAL(4)
      UE=0.0d0
      UC=XVAL(5)

      DA=XVAL(7)
      DB=XVAL(8)
      DE=0.0d0
      DC=XVAL(9)

      GA=XVAL(11)
      GB=XVAL(12)
      GE=0.0d0
      GC=XVAL(13)

      SN=XVAL(14)
      SA=XVAL(15)
      SB=XVAL(16)
      SE=0.0d0
      SC=XVAL(17)
      
      DLN=XVAL(18)
      DLA=XVAL(19)
c      DLB=XVAL(20)
      DLB=XVAL(16)+2.0d0
      DLE=0.0d0
      DLC=XVAL(21)
      AS=XVAL(1)
       CALL QNRSET('ALFAS',AS)
C--   Input quark distributions at Q2 = Q02 GeV2
C-- WRITE IN ACTUAL VALUES TO SAVE USING dgamma
C       UN=2./AREA(UA-1,UB,UE,UC)
C       DN=1./AREA(DA-1,DB,DE,DC)
      UN = XVAL(2)
      DN=XVAL(6)
c       UBDBN=DLN/AREA(DLA-1,DLB,DLE,DLC)
c      call grgive(nxgri,xxmin,xxmax,nqgri,qqmin,qqmax)
      nxgri=nx
      DO IX = 1,NXGRI
        xX = XFROMIX(IX)
        UPVAL=UN*FF_LHA(Xx,UA,UB,UE,UC)
	

        DNVAL=DN*FF_LHA(Xx,DA,DB,DE,DC)


        SEA=SN*FF_LHA(Xx,SA,SB,SE,SC)

C        GN=(1-UN*AREA(UA,UB,UE,UC)-
C     .        DN*AREA(DA,DB,DE,DC)-
C     .        SN*AREA(SA,SB,SE,SC))/AREA(GA,GB,GE,GC)
        GN=XVAL(10)
        GLUON=GN*FF_LHA(Xx,GA,GB,GE,GC)
	


c        UMD=UBDBN*FF_LHA(X,DLA,DLB,DLE,DLC)
        UMD=DLN*FF_LHA(Xx,DLA,DLB,DLE,DLC)
        SINGL=UPVAL+DNVAL+SEA        

c        print *,un,dn,sn,gn,dln

        CSEA=0.0
        SSEA=0.2*SEA
        USEA=(SEA-SSEA-CSEA)/2.-UMD
        DSEA=(SEA-SSEA-CSEA)/2.+UMD        
        UPLUS=UPVAL+USEA-1./NACT*SINGL
        DPLUS=DNVAL+DSEA-1./NACT*SINGL
        SPLUS=SSEA-1./NACT*SINGL


        CALL QNPSET('GLUON',IX,IQ0,GLUON)
        CALL QNPSET('SINGL',IX,IQ0,SINGL)
        CALL QNPSET('UPLUS',IX,IQ0,UPLUS)
        CALL QNPSET('DPLUS',IX,IQ0,DPLUS)
        CALL QNPSET('SPLUS',IX,IQ0,SPLUS)
        CALL QNPSET('UPVAL',IX,IQ0,UPVAL)
        CALL QNPSET('DNVAL',IX,IQ0,DNVAL)
      ENDDO
C--THINGS ARE FINE FOR HEAVY SO DO IT
      IF (HEAVY) THEN

C--     Evolve over whole grid

        CALL EVOLSG(IQ0,1,NQGRI)
        CALL EVPLUS('UPLUS',IQ0,1,NQGRI)
        CALL EVPLUS('DPLUS',IQ0,1,NQGRI)
        CALL EVPLUS('SPLUS',IQ0,1,NQGRI)
        CALL EVOLNM('UPVAL',IQ0,1,NQGRI)
        CALL EVOLNM('DNVAL',IQ0,1,NQGRI)

      ELSE

C--     Evolve gluon, singlet and valence over whole grid
      
        CALL EVOLSG(IQ0,1,NQGRI)
        CALL EVOLNM('UPVAL',IQ0,1,NQGRI)
        CALL EVOLNM('DNVAL',IQ0,1,NQGRI)
        
C--     Be more careful with plus distributions

        IF (NACT.EQ.3) THEN
C--THINGS ARE ALSO FINE IF 1Q0 IS BELOW 1QC THEN CLEARLY CSEA=0. IS OK
C--SITUATION CD BE 1<Q0<Q2C<Q2B ETC
C--GO DOWN QO TO 1 UP Q0 TO  Q2C
          CALL EVPLUS('UPLUS',IQ0,1,IQC)
          CALL EVPLUS('DPLUS',IQ0,1,IQC)
          CALL EVPLUS('SPLUS',IQ0,1,IQC)
C--DEAL WITH CHARM THRESH
          FACTOR=1./3.-1./4.
          CALL QADDSI('UPLUS',IQC,FACTOR)
          CALL QADDSI('DPLUS',IQC,FACTOR)
          CALL QADDSI('SPLUS',IQC,FACTOR)
          CALL QNPNUL('CPLUS')
          FACTOR=-1/4.
          CALL QADDSI('CPLUS',IQC,FACTOR)
C--GO UP Q2C TO Q2B
          CALL EVPLUS('UPLUS',IQC,IQC,IQB)
          CALL EVPLUS('DPLUS',IQC,IQC,IQB)
          CALL EVPLUS('SPLUS',IQC,IQC,IQB)
          CALL EVPLUS('CPLUS',IQC,IQC,IQB)

        ELSEIF(NACT.EQ.4)THEN
C--1<1QC<1Q0<1QB<ETC
C--NOW WE NEED A RETHINK OF THE INITIAL CONDITIONS
C--FIRST DEAL WITH CPLUS TRUNING ON AT IQC
C--AND GOING UP TO IQB (MIDDLE AGAIN)
          CALL QNPNUL('CPLUS')
          FACTOR=-1/4.
          CALL QADDSI('CPLUS',IQC,FACTOR)
          CALL EVPLUS('CPLUS',IQC,IQC,IQB)
C--REDIFINE THE PLUS DUSTNS TAKING INTO ACCOUNT CPLUS
          CALL QNPNUL('SPLUS')
          CALL QNPNUL('UPLUS')
          CALL QNPNUL('DPLUS')
C--NOW YOU NEED A DO LOOP AGAIN
          DO IX = 1,NXGRI
          Xx = XFROMIX(IX)
          UPVAL=UN*FF_LHA(Xx,UA,UB,UE,UC)

          DNVAL=DN*FF_LHA(Xx,DA,DB,DE,DC)
          SEA=SN*FF_LHA(Xx,SA,SB,SE,SC)
          SINGL=UPVAL+DNVAL+SEA
          UMD=DLN*FF_LHA(Xx,DLA,DLB,DLE,DLC)
          SSEA=0.2*SEA
          CSEA=QPDFIJ('CPLUS',IX,IQ0,IFL) + 1./NACT*SINGL
          USEA=(SEA-SSEA-CSEA)/2.-UMD
          DSEA=(SEA-SSEA-CSEA)/2.+UMD        
          UPLUS=UPVAL+USEA-1./NACT*SINGL
          DPLUS=DNVAL+DSEA-1./NACT*SINGL
          SPLUS=SSEA-1./NACT*SINGL
          CALL QNPSET('UPLUS',IX,IQ0,UPLUS)
          CALL QNPSET('DPLUS',IX,IQ0,DPLUS)
          CALL QNPSET('SPLUS',IX,IQ0,SPLUS)
          ENDDO
C-NOW DO MIDDLE FOR UPLUS,DPLUS,SPLUS
          CALL EVPLUS('UPLUS',IQ0,IQC,IQB)
          CALL EVPLUS('DPLUS',IQ0,IQC,IQB)
          CALL EVPLUS('SPLUS',IQ0,IQC,IQB)
C--THEN GO DOWN IQC TO 1
             IF(IQC.GT.1)THEN
             FACTOR=1./4.-1./3.
             CALL QADDSI('UPLUS',IQC,FACTOR)
             CALL QADDSI('DPLUS',IQC,FACTOR)
             CALL QADDSI('SPLUS',IQC,FACTOR)
             CALL EVPLUS('UPLUS',IQC,1,IQC)
             CALL EVPLUS('DPLUS',IQC,1,IQC)

             CALL EVPLUS('SPLUS',IQC,1,IQC)
             ENDIF

        ENDIF
C--THEN DEAL WITH B THRESHOLD FOR ALL4
        FACTOR=1./4.-1./5.
        CALL QADDSI('UPLUS',IQB,FACTOR)
        CALL QADDSI('DPLUS',IQB,FACTOR)
        CALL QADDSI('SPLUS',IQB,FACTOR)
        CALL QADDSI('CPLUS',IQB,FACTOR)
C--THEN GO UP
        IF(IQB.LT.NQGRI)THEN 
       CALL EVPLUS('UPLUS',IQB,IQB,NQGRI)
        CALL EVPLUS('DPLUS',IQB,IQB,NQGRI)
        CALL EVPLUS('SPLUS',IQB,IQB,NQGRI)
        CALL EVPLUS('CPLUS',IQB,IQB,NQGRI)
        ENDIF
cC--THEN DEAL WITH B TURNING ON AT IQB AND GO UP
        CALL QNPNUL('BPLUS')
        FACTOR=-1./5.
       CALL QADDSI('BPLUS',IQB,FACTOR)
        CALL EVPLUS('BPLUS',IQB,IQB,NQGRI)
      ENDIF  
c
      call getnset(iset)
      call save_pdfqcd(iset)

      return
*
      end
C------------------------------------------------------------------------ 
c      DOUBLE PRECISION FUNCTION AREA(A,B,E,C)
C------------------------------------------------------------------------ 
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c 
c      IF (A.LE.-0.99.OR.B.LE.-0.99) THEN
c        AREA=1D6
c        RETURN
c      ENDIF
c
c      AR1=(DGAMMA_LHA(A+1)*DGAMMA_LHA(B+1))/DGAMMA_LHA(A+B+2)
c      AR2=E*(DGAMMA_LHA(A+1.5)*DGAMMA_LHA(B+1))/DGAMMA_LHA(A+B+2.5)
c      AR3=C*(DGAMMA_LHA(A+2)*DGAMMA_LHA(B+1))/DGAMMA_LHA(A+B+3)
c      AREA=AR1+AR2+AR3
c     
c      IF (AREA.LE.1D-6) AREA=1D-6
c
c      RETURN
c      END
c
C------------------------------------------------------------------------ 
      DOUBLE PRECISION FUNCTION FF_LHA(X,A,B,E,C)
C------------------------------------------------------------------------ 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
      FF_LHA=X**A*(1D0-X)**B*(1+E*SQRT(X)+C*X)
     
      RETURN
      END
C------------------------------------------------------------------------ 
      SUBROUTINE DVZERO (A,N)
C
C CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71, modif. 24/05/87 to set integer zero
C                 modif. 25/05/94 to depend on QINTZERO
C
      implicit real*8 (a-h,o-z)
      DIMENSION A(*)
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= 0d0
      RETURN
      END
      



