#include "isajet/pilot.h"
      SUBROUTINE EEMAX
C          FIND UPPER BOUND FOR E+E- CROSS SECTION SUMMED OVER ALLOWED
C          TYPES.
C          VER 7.17: ENSURE XJMIN < XX < XJMAX
C          VER 7.42: ENACT BREMSSTRAHLUNG EFFECT
#if defined(CERNLIB_IMPNONE)
      IMPLICIT NONE
#endif
#include "isajet/itapes.inc"
#include "isajet/jetsig.inc"
#include "isajet/eepar.inc"
#include "isajet/jetlim.inc"
#include "isajet/jetpar.inc"
#include "isajet/primar.inc"
#include "isajet/xmssm.inc"
#include "isajet/sssm.inc"
#include "isajet/brembm.inc"
C
      REAL ETEST(3),XDI(3),RSH,XD,XDUMMY,SSFEL,DX,XX
      INTEGER NET,NXD,NX,NX1,IET,IXD,IX,I
C
      NET=1
      NXD=1
      ETEST(1)=ECM
C     Initialize beam/brem spectra convolution and fit
      IF (IBEAM) THEN
        QSQBM=RSHMAX**2
        EB=RSHMAX/2.
        WRITE(ITLIS,*) '   BEGINNING BREM/BEAM CONVOLUTION AND FIT...'
        XDUMMY=SSFEL(.1,1)
      END IF
      IF (IBREM) THEN
        NET=3
        NXD=3
        ETEST(1)=RSHMIN
        IF (RSHMAX.GT.AMZ.AND.RSHMIN.LT.AMZ) THEN
          ETEST(2)=AMZ
        ELSE
          ETEST(2)=RSHMIN+(RSHMAX-RSHMIN)/2.
        END IF
        ETEST(3)=MIN(RSHMAX,.999*ECM)
      END IF
      SGMXEE=0.
      NX=50
      IF(FIXYJ(1)) NX=1
      NX1=NX+1
      DX=(XJMAX(1)-XJMIN(1))/NX
C          SCAN IN X=COS(THETA)
      DO IET=1,NET
       RSH=ETEST(IET)
       SHAT=RSH*RSH
       XDI(1)=-(1.-SHAT/SCM)
       XDI(2)=0.
       XDI(3)=-XDI(1)
       DO IXD=1,NXD
        XD=XDI(IXD)
        X1=(XD+SQRT(XD**2+4*SHAT/SCM))/2.
        X2=X1-XD
        DO 100 IX=1,NX1
         XX=XJMIN(1)+DX*(IX-1)
         IF(XX.LT.XJMIN(1)) XX=XJMIN(1)
         IF(XX.GT.XJMAX(1)) XX=XJMAX(1)
         CTH(1)=XX
         CTH(2)=-XX
         DO 110 I=1,2
         XJ(I)=CTH(I)
         TH(I)=ACOS(CTH(I))
         STH(I)=SIN(TH(I))
         PT(I)=HALFE*STH(I)
         IF(IX.EQ.1) YJ(I)=YJMIN(I)
         IF(IX.EQ.NX1) YJ(I)=YJMAX(I)
         IF(IX.GT.1.AND.IX.LT.NX1)
     1   YJ(I)=.5*ALOG((1.+CTH(I))/(1.-CTH(I)))
110      CONTINUE
C          COMPUTE CROSS SECTION
         IF (GOMSSM) THEN
           CALL SIGSSE
         ELSE
           CALL SIGEE
         END IF
         IF(SIGMA.GT.SGMXEE) SGMXEE=SIGMA
100     CONTINUE
       END DO
      END DO
C          REQUIRE CROSS SECTION BE POSITIVE
      WRITE(ITLIS,1000) SGMXEE
1000  FORMAT(///' MAXIMUM D(SIGMA)/D(COS THETA) = ',E12.4)
      IF(SGMXEE.GT.0) RETURN
      STOP 99
      END
