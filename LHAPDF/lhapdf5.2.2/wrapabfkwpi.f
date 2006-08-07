      subroutine ABFKWPevolve(xin,qin,pdf)
      include 'parmsetup.inc'
      PARAMETER(NX=50)
      PARAMETER(NQ=19)
      real*8 xin,qin,pdf(-6:6),xval(45),qcdl4,qcdl5  
      real*8 upv,dnv,usea,dsea,str,chm,bot,top,glu
      real*8 calcpi(8,20,25,3),calcpio(8,20,25),parpi(40,3)
      common /ABFKWP/ CALCPI,CALCPIO,PARPI,lastmem 
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mmem
      common/NAME/name,nmem,ndef,mmem
      integer nset
      save 
      
      iimem = imem
      if(iimem.eq.0) iimem = 1
      if(iimem.le.3) then
       call ABFKWxx(iimem,xin,qin,upv,dnv,usea,dsea, str,chm,glu)
      endif  
      

      pdf(-6)= 0.0d0
      pdf(6)= 0.0d0
      pdf(-5)= 0.0d0
      pdf(5 )= 0.0d0
      pdf(-4)= chm
      pdf(4 )= chm
      pdf(-3)= str
      pdf(3 )= str
      pdf(-2)= usea
      pdf(2 )= upv+usea
      pdf(-1)= dsea
      pdf(1 )= dnv+dsea
      pdf(0 )= glu
      
      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry ABFKWPread(nset)
      read(1,*)nmem(nset),ndef(nset)
c      print *,nmem,ndef
      lastmem = -999
      do j=1,3
        read(1,*)(parpi(k,j),k=1,4)
        read(1,*)(parpi(k,j),k=5,8)
        read(1,*)(parpi(k,j),k=9,12)
        read(1,*)(parpi(k,j),k=13,16)
        read(1,*)(parpi(k,j),k=17,20)
        read(1,*)(parpi(k,j),k=21,24)
        read(1,*)(parpi(k,j),k=25,28)
        read(1,*)(parpi(k,j),k=29,32)
        read(1,*)(parpi(k,j),k=33,36)
        read(1,*)(parpi(k,j),k=37,40)
        do l=1,25
          do k=1,20
            read(1,*)(CALCPI(m,k,l,j),m=1,4)
            read(1,*)(CALCPI(m,k,l,j),m=5,8)
          enddo
        enddo
      enddo
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry ABFKWPalfa(alfas,qalfa)
        call getnset(iset)
	call GetOrderAsM(iset,iord)
        call Getlam4M(iset,imem,qcdl4)
        call Getlam5M(iset,imem,qcdl5)
        call aspdflib(alfas,Qalfa,iord,qcdl5)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry ABFKWPinit(Eorder,Q2fit)
      return
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry ABFKWPpdf(mem)
      imem = mem
      return
c
 1000 format(5e13.5)
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:27:26  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:28:53  plothow
* Version 7.01
*
*
      SUBROUTINE ABFKWxx(imem,DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL)
      double precision
     +       PARPI(40,3),CALCPI(8,20,25,3),CALCPIO(8,20,25),ZEROD,
     +       DX,DQ,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DGL
      REAL    X, Q, UPV, DNV, USEA, DSEA, STR, CHM, GL
 
      common /ABFKWP/CALCPI,CALCPIO,PARPI,lastmem 
     
c      COMMON/W5051Ixx/CALCPIO
      REAL   XPDF(7)
      DATA ZEROD/0.D0/
C----------------------------------------------------------------------
       DATA ISTART/0/
       SAVE ISTART,OWLAM2,Q02PI
C
      if(imem.ne.lastmem) then
         istart = 0
	 lastmem = imem
      endif
      IF (ISTART.EQ.0) THEN
        ISTART=1
        DO 11 K=1,25
        DO 11 I=1,20
        DO 11 M=1,8
   11   CALCPIO(M,I,K) = CALCPI(M,I,K,imem)
           OWLAM=PARPI(1,imem)
           OWLAM2=OWLAM**2
           Q02PI=PARPI(39,imem)
           Q2MAX=PARPI(40,imem)
         ENDIF
C
C the conventions are : q(1)=x*u, q(2)=x*d, q(3)=x*str, q(4)=x*usea,
C                       q(5)=x*dsea, q(6)=x*charm, q(7)=x*gluon
C
      X = DX
      Q = DQ
      Q2 = Q*Q
      IDQ2=2
      SB=0.
      IF(Q2-Q02PI) 1,1,2
    2 IF(IDQ2-1) 1,1,3
    3 SB= LOG( LOG( MAX(Q02PI,Q2)/OWLAM2)/ LOG(Q02PI/OWLAM2))
    1 CALL AURPIx(1,0,X,SB,XPDF(1))
      CALL AURPIx(2,0,X,SB,XPDF(2))
      CALL AURPIx(3,0,X,SB,XPDF(3))
      CALL AURPIx(4,0,X,SB,XPDF(4))
      CALL AURPIx(5,0,X,SB,XPDF(5))
      CALL AURPIx(8,0,X,SB,XPDF(6))
      CALL AURPIx(7,0,X,SB,XPDF(7))
C
      DUPV=XPDF(1) - XPDF(4)
      DDNV=XPDF(2) - XPDF(5)
      DUSEA=XPDF(4)
      DDSEA=XPDF(5)
      DSTR=XPDF(3)
      DCHM=XPDF(6)
      DGL =XPDF(7)
C
      RETURN
      END
c==============================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:27:36  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:03  plothow
* Version 7.01
*
*
C
      SUBROUTINE AURPIx(I,NDRV,X,S,ANS)
      double precision
     +       CALCPI(8,20,25,3),CALCPIO(8,20,25),parpi(40,3)
      common /ABFKWP/CALCPI,CALCPIO,parpi,lastmem
c      COMMON/W5051I4/CALCPIO
      REAL   F1(25),F2(25)
      DATA DELTA/.10/
      ANS=0.
      IF(X.GT.0.9985) RETURN
      IF(I.EQ.3.AND.X.GT.0.95) RETURN
      IF(I.EQ.8.AND.X.GT.0.95) RETURN
      IS=S/DELTA+1
      IS1=IS+1
      DO 1 L=1,25
      KL=L+NDRV*25
      F1(L)=CALCPIO(I,IS,KL)
      F2(L)=CALCPIO(I,IS1,KL)
    1 CONTINUE
      A1=AUGETFV(X,F1)
      A2=AUGETFV(X,F2)
      S1=(IS-1)*DELTA
      S2=S1+DELTA
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)
      RETURN
      END
c===============================================================
*
* $Id$
*
* $Log$
* Revision 1.2  2005/10/07 15:15:05  whalley
* Changes to most files for V5 - multiset initializations
*
* Revision 1.1.1.1  2005/05/06 14:54:43  whalley
* Initial CVS import of the LHAPDF code and data sets
*
* Revision 1.1.1.2  1996/10/30 08:27:34  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:02  plothow
* Version 7.01
*
*
C
      FUNCTION AUGETFV(X,FVL)
C  LOGARITHMIC INTERPOLATOR - WATCH OUT FOR NEGATIVE
C  FUNCTIONS AND/OR X VALUES OUTSIDE THE RANGE 0 TO 1.
C  NOTE: DIMENSION OF FVL IS OVERWRITTEN BY VALUE USED
C  IN MAIN ROUTINE.
      DIMENSION FVL(25),XGRID(25)
      DATA NX,XGRID/25,.001,.002,.004,.008,.016,.032,.064,.1,.15,
     *.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95/
      AUGETFV=0.
      DO 1 I=1,NX
      IF(X.LT.XGRID(I)) GO TO 2
    1 CONTINUE
    2 I=I-1
      IF(I.EQ.0) THEN
         I=I+1
      ELSE IF(I.GT.23) THEN
         I=23
      ENDIF
      J=I+1
      K=J+1
      AXI= LOG(XGRID(I))
      BXI= LOG(1.-XGRID(I))
      AXJ= LOG(XGRID(J))
      BXJ= LOG(1.-XGRID(J))
      AXK= LOG(XGRID(K))
      BXK= LOG(1.-XGRID(K))
      FI= LOG(ABS(FVL(I)) +1.E-15)
      FJ= LOG(ABS(FVL(J)) +1.E-16)
      FK= LOG(ABS(FVL(K)) +1.E-17)
      DET=AXI*(BXJ-BXK)+AXJ*(BXK-BXI)+AXK*(BXI-BXJ)
      ALOGA=(FI*(AXJ*BXK-AXK*BXJ)+FJ*(AXK*BXI-AXI*BXK)+FK*(AXI*BXJ-AXJ*
     $ BXI))/DET
      ALPHA=(FI*(BXJ-BXK)+FJ*(BXK-BXI)+FK*(BXI-BXJ))/DET
      BETA=(FI*(AXK-AXJ)+FJ*(AXI-AXK)+FK*(AXJ-AXI))/DET
      IF(ABS(ALPHA).GT.99..OR.ABS(BETA).GT.99..OR.ABS(ALOGA).GT.99.)
     1RETURN
C      IF(ALPHA.GT.50..OR.BETA.GT.50.) THEN
C         WRITE(6,2001) X,FVL
C 2001    FORMAT(8E12.4)
C         WRITE(6,2001) ALPHA,BETA,ALOGA,DET
C      ENDIF
      AUGETFV=EXP(ALOGA)*X**ALPHA*(1.-X)**BETA
      RETURN
      END
c============================================================
