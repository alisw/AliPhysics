*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:32:51  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:17  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM TESTGAMO
C
#include "pdf/impdp.inc"
      PARAMETER (LFLPRT = 2)
C+SEQ, W5051P2.
      PARAMETER (MODEMX = 281)
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
C
#include "pdf/expdp.inc"
     +       DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL
      REAL       GLUCK
      COMMON/W5051I8/GLUCK
#include "pdf/expdp.inc"
     +               XI(32),WI(32),XX(33)
      COMMON/W5051I9/XI,WI,XX,NTERMS
#include "pdf/expdp.inc"
     +       XV(3),ZQ(5)
      DATA XV/0.05D0,0.50D0,0.95D0/
      DATA ZQ/10.D0,100.D0,1000.D0,1.D+4,5.5D+5/
      DATA FLAV/4.D0/
      DATA ALQED/137.D0/
      DATA ICALL/0/
      CHARACTER*20 PARM(20)
#include "pdf/expdp.inc"
     +       VAL(20)
C
      N6 = L6
      IFLPRT = LFLPRT
C
      IF(GLUCK.LT.0.01) GLUCK = 0.
      PI=4.D0*ATAN(1.D0)
      PI2=PI**2
      CF=4.D0/3.D0
      NF=INT(FLAV+0.1D0)
      FLAVOR=FLAV
      B0=11.D0-2.D0*NF/3.D0
      B1=102.D0-38.D0*NF/3.D0
      CCOEG=2.D0/9.D0
      IF(NF.EQ.4) CCOEG=5.D0/18.D0
      COEG=2.D0*FLAVOR*CCOEG
      CALL WATE32
C
      DO 300 ISTRUC=272,272
      ICALL = 0
      FIRST=.TRUE.
      VAL(1)  = ISTRUC
      PARM(1) = 'Mode'
      IF(ISTRUC.GE.201 .AND. ISTRUC.LE.205) GOTO 300
      IF(ISTRUC.GE.208 .AND. ISTRUC.LE.209) GOTO 300
      IF(ISTRUC.GE.214 .AND. ISTRUC.LE.219) GOTO 300
      IF(ISTRUC.GE.223 .AND. ISTRUC.LE.229) GOTO 300
      IF(ISTRUC.GE.233 .AND. ISTRUC.LE.270) GOTO 300
      IF(ISTRUC.GE.274 .AND. ISTRUC.LE.279) GOTO 300
      IF(ISTRUC.GE.282                    ) GOTO 300
C
         CALL PDFSET(PARM,VAL)
         Q02 = Q2MIN
         OWLAM = QCDL4
         OWLAM2=OWLAM**2
         IORD = LO - 1
*
      DO 100 J=3,3
      Q2=Q02
      IF(J.GE.2) Q2=ZQ(J-1)
      ALQ2= LOG(Q2/OWLAM2)
      ALFPI= 2.D0 /(B0*ALQ2+B1* LOG(ALQ2)/B0)
      DO 70 K=1,3
      X=XV(K)
      DX = X
      Q = SQRT(Q2)
      DQ = Q
      CALL STRUCTM(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
      UB =DUB
      DB =DDB
      SB =DSB
      CB =DCB
      GL =DGL
      F2=(4.D0*  DUB   +  DDB   +  DSB)/9.D0*2.D0
      F2C=F2+(4.D0/9.D0)*DCB*2.D0
      IF(ICALL.EQ.0) THEN
        ICALL = 1
        WRITE (N6,1000) Q02,IORD,OWLAM,GLUCK,FLAV
 1000   FORMAT(1X,' Q02= ',F7.3,'   IORD= ',I2,'   LAMBDA=',F7.4,
     +               '   GLUCK= ',F5.2,'   NF = ',F5.2)
        ENDIF
      IF (K.EQ.1) THEN
        WRITE (N6,1001) Q2
 1001   FORMAT(/,'  Q2=  ',E10.3,' GEV**2')
        WRITE(N6,1002)
 1002   FORMAT(/,4X,' X',6X,'X*UB  ',6X,'X*DB  ',6X,'X*SB  ',6X,'X*CB  '
     +                  ,6X,'X*GL',8X,' F2 ',6X,' F2/ALPHA'
     +                            ,5X,' F2C ',5X,' F2C/ALPHA')
      ENDIF
      IF (GLUCK.LT..5) GOTO 5
      F2=DUV
 5    IF(IORD)13,13,14
14    CONTINUE
      AL=0.5D0*ALFPI
      AL1= LOG(1.D0-X)
      XQ=F2+F2*AL*CF*(-9.D0-2.D0*PI2/3.D0+AL1*(-3.D0+2.D0*AL1))
C     WRITE(N6,1111) X,Q,F2,AL1,XQ
C1111 FORMAT(1X,' X,F2 = ',5E13.5)
      DO 81 I=1,NTERMS
      Y=0.5D0*(1.D0-X)*XI(I)+0.5D0*(1.D0+X)
      XY=X/Y
      AL1= LOG(1.D0-Y)
      DX = XY
      CALL STRUCTM(DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
      XQQ=(4.D0*  DUB   +  DDB   +  DSB)/9.D0*2.D0
      IF (GLUCK.LT..5) GOTO 6
      XQQ=DUV
 6    C22=CF*(6.D0+4.D0*Y-2.D0*(1.D0+Y*Y)/
     +          (1.D0-Y)* LOG(Y)-2.D0*(1.D0+Y)*AL1)
      C23=CF*(-3.D0+4.D0*AL1)/(1.D0-Y)
      CG2=(-1.D0+8.D0*Y*(1.D0-Y)+(1.D0-2.D0*Y+2.D0*Y*Y)*
     +                   LOG(1.D0/Y-1.D0))*COEG
      XQ=XQ +0.5D0*(1.D0-X)*WI(I)*AL*(C22*XQQ+C23*(XQQ- F2))
      XQ=XQ +0.5D0*(1.D0-X)*WI(I)*AL*CG2*DGL
C     WRITE(N6,1112) Y,AL1,C22,C23,CG2,XQQ,XQ
C1112 FORMAT(1X,' Y,C22 = ',7E13.5)
 81   CONTINUE
      F2=3.D0*3.D0*(2.D0/27.D0)*(1.D0/(2.D0*PI*ALQED))*2.D0
      F2=XQ+F2*(8.D0*X*(1.D0-X)-1.D0+(1.D0-2.D0*X+2.D0*X*X)*
     +            LOG(1.D0/X-1.D0))*X
      F2C=F2+(4.D0/9.D0)*CB*2.D0
   13 F2N=F2*ALQED
      F2CN=F2C*ALQED
      WRITE(N6,1003)X,UB,DB,SB,CB,GL,F2,F2N,F2C,F2CN
 1003 FORMAT(1X,' ',F6.3,9E12.4)
   70 CONTINUE
  100 CONTINUE
      WRITE(N6,*) '  '
  300 CONTINUE
      STOP
      END
