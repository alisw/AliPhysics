*
* $Id$
*
* $Log$
* Revision 1.1.1.4  2000/04/19 13:30:01  mclareni
* Import version 8.04
*
* Revision 1.4  2000/04/09 10:06:56  plothow
* new version 8.00 - MRS 98 and 99, CTEQ5 and GRV98 added
*
* Revision 1.3  1997/07/02 15:03:11  pdflib
* potential bug removed in GRSGALO
*
* Revision 1.2  1996/10/28 16:10:10  pdflib
* PDFLIB new version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:17  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      PROGRAM TESTGAM
C
#include "pdf/impdp.inc"
      PARAMETER (LFLPRT = 2)
#include "pdf/w5051p2.inc"
#include "pdf/w5051p7.inc"
#include "pdf/w50510.inc"
#include "pdf/w50511.inc"
#include "pdf/w50512.inc"
#include "pdf/w50513.inc"
#include "pdf/w50516.inc"
#include "pdf/w50517.inc"
#include "pdf/w505110.inc"
#include "pdf/w505120.inc"
#include "pdf/w505122.inc"
C
#include "pdf/expdp.inc"
     +       DX,DQ,DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL,
     +       DQ2,DP2
#include "pdf/expdp.inc"
     +               XI(32),WI(32),XX(33)
      COMMON/W5051I9/XI,WI,XX,NTERMS
#include "pdf/expdp.inc"
     +       XV(3),ZQ(5)
      CHARACTER*20 PARM(NCHDIM)
#include "pdf/expdp.inc"
     +       VAL(NCHDIM)
C
      DATA XV/0.05D0,0.50D0,0.95D0/
      DATA ZQ/10.D0,100.D0,1000.D0,1.D+4,5.5D+5/
      DATA P2/5./
      DATA IP2/0/
      DATA FLAV/4.D0/
      DATA ALQED/137.D0/
      DATA ICALL/0/
C
      N6 = L6
      IFLPRT = LFLPRT
C
C*****adjustment of the VDM contribution*******************************
      KA=1.
C*****mass of the charm quark******************************************
      CM=1.41
C*****CHARME=0. -> no charm contribution ******************************
      CHARME=1.
C
      PI=4.D0*ATAN(1.D0)
      PI2=PI**2
      CF=4.D0/3.D0
      NF=INT(FLAV+0.1D0)
      IF(NF.EQ.0) NF=3
      FLAVOR=FLAV
      B0=11.D0-2.D0*NF/3.D0
      B1=102.D0-38.D0*NF/3.D0
      CCOEG=2.D0/9.D0
C      IF(NF.EQ.4) CCOEG=5.D0/18.D0
C      COEG=2.D0*FLAVOR*CCOEG
      COEG=2.D0*FLAVOR*CCOEG
      CALL WATE32
C
C First call to PDFSET to initialize COMMON/W505120/
      PARM(1) = 'Init0'
      VAL(1)  = 0.D0
      CALL PDFSET(PARM,VAL)
C.
C Test new PDFLIB versions (4.0 or bigger) for Photons only
      NPTYPE = 3
         DO KGRO = 1,NGRMAX
            IF(NPGSMX(NPTYPE,KGRO) .EQ. 0) GOTO 500
            DO ISET=1,NPGSMX(NPTYPE,KGRO)
               WRITE (N6,*) ' ==========================================
     +============================================'
               FIRST=.TRUE.
               PARM(1) = 'Nptype'
               VAL(1)  = NPTYPE
               PARM(2) = 'Ngroup'
               VAL(2)  = KGRO
               PARM(3) = 'Nset'
               VAL(3)  = ISET
               CALL PDFSET(PARM,VAL)
C
         ICALL = 0
         Q02 = Q2MIN
         OWLAM = QCDL4
         OWLAM2=OWLAM**2
         IORD = LO - 1
*
      DO 100 J=3,3
      Q2=Q02
      IF(J.GE.2) Q2=ZQ(J-1)
      XTH=1.D0/(1.D0+4.D0*CM**2/Q2)
      ALQ2= LOG(Q2/OWLAM2)
      ALFPI= 2.D0 /(B0*ALQ2+B1* LOG(ALQ2)/B0)
      DO 70 K=1,3
      X=XV(K)
      DX = X
      Q = SQRT(Q2)
      DQ = Q
      DQ2 = Q2
      DP2 = P2
      IP = IP2
      CUT=1.
      IF(CHARME.LT.0.5) CUT=0.
      CUTG=1.
      IF(X.GT.XTH) CUT=0.
      CALL STRUCTP(DX,DQ2,DP2,IP,
     +             DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
      UB =DUB
      DB =DDB
      SB =DSB
      CB =DCB
      GL =DGL
      F2=(4.D0*  DUB   +  DDB   +  DSB)/9.D0
      F2C=4.D0/9.D0*(FLAVOR-3.D0)*DCB*CUT
      F2=F2+F2C
      IF(ICALL.EQ.0) THEN
        ICALL = 1
        WRITE (N6,1000) Q02,P2,CHARME,CUT,IORD,OWLAM,FLAV
 1000   FORMAT(1X,' Q02= ',F7.3,'   P2= ',F7.3,
     +               '   CHARM= ',F3.1,'   CUT= ',F3.1,
     +               '   IORD= ',I2,'   LAMBDA=',F7.4,'   NF = ',F5.2)
        ENDIF
      IF (K.EQ.1) THEN
        WRITE (N6,1001) Q2
 1001   FORMAT(/,'  Q2=  ',E10.3,' GEV**2')
        WRITE(N6,1002)
 1002   FORMAT(/,4X,' X',6X,'X*UB  ',6X,'X*DB  ',6X,'X*SB  ',6X,'X*CB  '
     +                  ,6X,'X*GL',8X,' F2 ',6X,' F2/ALPHA'
     +                            ,5X,' F2C ',5X,' F2C/ALPHA')
      ENDIF
      IF(IORD)13,13,14
14    CONTINUE
      AL=0.5D0*ALFPI
      AL1= LOG(1.D0-X)
      XQ=F2+F2*AL*CF*(-9.D0-2.D0*PI2/3.D0+AL1*(-3.D0+2.D0*AL1))
C     WRITE(N6,1111) X,Q,F2,AL1,XQ
C1111 FORMAT(1X,' X,F2 = ',5E13.5)
      XM=1.-X
      XP=1.+X
      DIS2=(XM**2-4.*X*XM*CM**2/Q2)*CUT
      DIS=SQRT(DIS2)
      YP=(XP+DIS)/2.
      YM=(XP-DIS)/2. 
      CUTY=1. 
      DO 81 I=1,NTERMS
      Y=0.5D0*(1.D0-X)*XI(I)+0.5D0*(1.D0+X)
      XY=X/Y
      AL1= LOG(1.D0-Y)
      DX = XY
      CALL STRUCTP(DX,DQ2,DP2,IP,
     +             DUV,DDV,DUB,DDB,DSB,DCB,DBB,DTB,DGL)
      XQQ=(4.D0*  DUB   +  DDB   +  DSB)/9.D0
      XQQ=XQQ+4.D0/9.D0*(FLAVOR-3.D0)*DCB*CUT*CUTY
      F2A=F2+F2C*(CUTY-1.)
 6    C22=CF*(6.D0+4.D0*Y-2.D0*(1.D0+Y*Y)/
     +          (1.D0-Y)* LOG(Y)-2.D0*(1.D0+Y)*AL1)
      C23=CF*(-3.D0+4.D0*AL1)/(1.D0-Y)
      COEGM=COEG+(8./9.)*(FLAVOR-3.)*CUTG*CUT
      CG2=(-1.D0+8.D0*Y*(1.D0-Y)+(1.D0-2.D0*Y+2.D0*Y*Y)*
     +                   LOG(1.D0/Y-1.D0))*COEGM
      XQ=XQ +0.5D0*(1.D0-X)*WI(I)*AL*(C22*XQQ+C23*(XQQ- F2A))
      XQ=XQ +0.5D0*(1.D0-X)*WI(I)*AL*CG2*(DGL+KA*GL)
C     WRITE(N6,1112) Y,AL1,C22,C23,CG2,XQQ,XQ
C1112 FORMAT(1X,' Y,C22 = ',7E13.5)
      CUTG=1.
   81 CUTY=1.
      COE3=3.D0*3.D0*(2.D0/27.D0)*(1.D0/(2.D0*PI*ALQED))*2.D0
      COE4=3.D0*(34.D0/81.D0)*(1.D0/(2.D0*PI*ALQED))*2.D0
      IF(NF.EQ.3) COE4=COE3
         F2=XQ+COE4*(8.D0*X*(1.D0-X)-1.D0+(1.D0-2.D0*X+2.D0*X*X)*
     +            LOG(1.D0/X-1.D0))*X
         F2M=XQ+COE3*(8.D0*X*(1.D0-X)-1.D0+(1.D0-2.D0*X+2.D0*X*X)*
     +            LOG(1.D0/X-1.D0))*X
      CMS=CM**2
      SC=Q2*(1.D0-X)/X
      BE=4.*CMS/SC
      IF(BE.GE.1.) WCM=0.
      IF(BE.GE.1.) GO TO 91
      SQ=SQRT(1.-BE)
      A1=((1.D0+SQ)/2.D0)**2 
      A2=(1.D0-X)/X
      WCM=(8.D0*(1.D0-X)*X-1.D0)*SQ+(X**2+(1.D0-X)**2)*LOG(A1*A2)
      WCM=3.D0*(4.D0/9.D0)/(PI*ALQED)*X*WCM
91    CONTINUE
      F2M=F2M+(4.D0/9.D0)*(FLAVOR-3.D0)*WCM*CUT
         F2C=F2M
   13    F2N=F2*ALQED
         F2CN=F2C*ALQED
         WRITE(N6,1003)X,UB,DB,SB,CB,GL,F2,F2N,F2C,F2CN
 1003    FORMAT(1X,' ',F6.3,9E12.4)
   70    CONTINUE
  100    CONTINUE
       ENDDO
       WRITE(N6,*) '  '
  500  CONTINUE
      ENDDO
      STOP
      END
