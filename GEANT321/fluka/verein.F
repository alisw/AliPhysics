*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE VEREIN.FOR
*COPY VEREIN
*
*=== verein ===========================================================*
*
      SUBROUTINE VEREIN(IT,LA,LT,RER,REL,RPXR,RPYR,RPZR,RPXL,RPYL,RPZL,
     *KR1R,KR2R,KR1L,KR2L,IHAD,LL,KFR1,KFR2,IMPS,IMVE,IB08,IA08,
     *IB10,IA10,B3,AS,B8,IAR,KFA1,KFA2,KFA3,KFA4,IOPT)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*  Verein89: slight revision by A. Ferrari                             *
*----------------------------------------------------------------------*
*
#include "geant321/finpar2.inc"
#include "geant321/part.inc"
      DIMENSION KFR1(*),KFR2(*),IMPS(6,6)
*
      PARAMETER (KMXJCM = 100)
      DIMENSION IV(KMXJCM)
      DIMENSION IMVE(6,6),IB08(6,21),IA08(6,21),IB10(6,21),IA10(6,21)
      REAL RNDM(1)
C*****VEREIN COMBINES THE TWO JETS INTO A COMPLETE E+E- EVENT OR INTO
C*****A COMPLETE QQQ (AQ,AQ,AQ) EVENT
C*****RER,REL,RPXR,RPXL,RPYR,RPYL,RPZR,RPZL ARE REST JET ENERGIES AND
C*****MOMENTA , KR1R,KR1L,KR2R,KR2L ARE REST JET FLAVOURS
*
*     R is for right, L for left!
*
*  Initialize Iv
*
      DATA IV /KMXJCM*0/
      IF(LT.EQ.0) GO TO 300
      WRITE(LUNOUT,202)IT
  202 FORMAT(1H0,I5,3H=IT)
  300 CONTINUE
C*****ONLY ONE RESONANCE WILL BE CREATED IF IT=0
      IF(IT.EQ.0) GO TO 100
      IF(KR2R.EQ.0.OR.KR1R.EQ.0)GO TO 10
      IF(KR1L.EQ.0.OR.KR2L.EQ.0)GO TO 10
C*****RESTJET CONTAINS FOUR QUARKS,THEREFORE TWO MESONS WILL BE FORMED
*
*  Reg is et equal to the sum of the left and right jet energy rest,
*  the same for momenta
*
      REG=RER+REL
      RPXG=RPXR+RPXL
      RPYG=RPYR+RPYL
      RPZG=RPZR+RPZL
      J=IT+1
      I=IT+2
      IT=J
*
*  Hklass uses always only Iv(i) with i >= it, so from the previous
*  instructions, it is always > than the original it , so the used
*  Iv's should be always 0 also in the bamjet original array ?????
*
      CALL HKLASS(IT,LT,LA,LL,KFR1,KFR2,KR1R,KR2R,KR1L,KR2L,IV,IMPS,
     &IMVE,IB08,IA08,IB10,IA10,AS,B8,KFA1,KFA2,KFA3,KFA4,IOPT)
      IHAD=IHAD+2
      IF(LT.EQ.1) WRITE(LUNOUT,301)IHAD
  301 FORMAT(1H0,13HHADRONANZAHL=,I5)
      GO TO 20
C*****RESTJET CONTAINS TWO OR THREE QUARKS,THEREFORE ONE MESON OR ONE
C*****BARYON OR ONE ANTIBARYON  WILL BE FORMED
   10 CONTINUE
      IAK=IT
      CALL GRNDM(RNDM,1)
      X=RNDM(1)
      IF(X.LE.0.5D0.AND.IAR.NE.0) IAK=IAR
      REG=RER+REL+HEF(IAK)
      RPXG=RPXR+RPXL+PXF(IAK)
      RPYG=RPYR+RPYL+PYF(IAK)
      RPZG=RPZR+RPZL+PZF(IAK)
      I=IT+1
      J=IAK
      IT=I
*
*  Hklass uses always only Iv(i) with i >= it, so from the previous
*  instructions, it is always > than the original it , so the used
*  Iv's should be always 0 also in the bamjet original array ?????
*
      CALL HKLASS(IT,LT,LA,LL,KFR1,KFR2,KR1R,KR2R,KR1L,KR2L,IV,IMPS,
     &IMVE,IB08,IA08,IB10,IA10,AS,B8,KFA1,KFA2,KFA3,KFA4,IOPT)
*  Hp1 is not used, also Hp2 !!!
*     HP2=HP1
      IHAD=IHAD+1
      IF(LT.EQ.1) WRITE(LUNOUT,301)IHAD
   20 CONTINUE
      RPG=SQRT(RPXG**2+RPYG**2+RPZG**2)
      IF(RPG.GT.REG) LA=3
      IF(LA.EQ.3.AND.LT.GT.0) WRITE(LUNOUT,71)
   71 FORMAT(2X,'REST JET MOMENTUM IS GREATER THAN REST JET ENERGY')
      IF(LA.EQ.3) RETURN
      RMG=SQRT((REG-RPG)*(REG+RPG))
      HM1=AMF(J)
      HM2=AMF(I)
      IF(RMG.LT.HM1+HM2) LA=3
      IF(LA.EQ.3.AND.LT.GT.0) WRITE(LUNOUT,71)
      IF(LA.EQ.3) RETURN
      IF(LT.EQ.0) GO TO 30
      WRITE(LUNOUT,60)REG,RMG,RPG,HM1,HM2,J,I
   60 FORMAT(1H0,19HEG,MG,PG,M1,M2,J,I=,5F8.4,2I3)
C     LORENZ PARAMETERS
   30 CONTINUE
      GAA=REG/RMG
      GABE=RPG/RMG
C     DECAY INTO 2 HADRONS IN THE CMS OF THE REMAINDER
      HE1=(RMG**2+HM1**2-HM2**2)/(2.D0*RMG)
      HE2=RMG-HE1
*   Hp1 is not used !
*     HP1=SQRT(HE1**2-HM1**2)
C     SAMPLES THE MOMENTUM DIRECTIONS OF THE LAST PARTICLES
      HE=HE1
      HMA=HM1
      CALL FKIMPU(HE,HMA,HPS,HPX,HPY,HPZ,LT,LL,B3)
      HP1X=HPX
      HP1Y=HPY
      HP1Z=HPZ
      HP2X=-HP1X
      HP2Y=-HP1Y
      IF (IOPT.EQ.4.AND.HM2.GT.HM1)  GO TO 999
      HP1Z=-HPZ
      HP2Z=HPZ
  999 CONTINUE
      HP2Z=-HP1Z
      IF(RPG.EQ.0.D0)GO TO 80
C     ROTATION BACK TO THE CMS*
      HEX=HE1*GAA+HP1Z*GABE
      HEY=HE2*GAA+HP2Z*GABE
      HP1Z=HE1*GABE+HP1Z*GAA
      HP2Z=HE2*GABE+HP2Z*GAA
      HE1=HEX
      HE2=HEY
C     ROTATION INTO THE CMS OF THE E+ E-  COLLISION
      X=HP1X
      Y=HP1Y
      Z=HP1Z
      COTE=RPZG/RPG
      SITE=SQRT((1.D0-COTE)*(1.D0+COTE))
      IF(SITE.EQ.0.D0) GO TO 11
      COV=RPXG/(RPG*SITE)
      SIV=RPYG/(RPG*SITE)
      GO TO 12
   11 CONTINUE
      SIV=1.D0
      COV=0.D0
   12 CONTINUE
      COPS=-SIV
      SIPS=COV
      CALL DRELAB(X,Y,Z,COTE,SITE,COPS,SIPS)
      HP1X=X
      HP1Y=Y
      HP1Z=Z
      X=HP2X
      Y=HP2Y
      Z=HP2Z
      CALL DRELAB(X,Y,Z,COTE,SITE,COPS,SIPS)
      HP2X=X
      HP2Y=Y
      HP2Z=Z
   80 CONTINUE
      HEF(J)=HE1
      PXF(J)=HP1X
      PYF(J)=HP1Y
      PZF(J)=HP1Z
      HEF(I)=HE2
      PXF(I)=HP2X
      PYF(I)=HP2Y
      PZF(I)=HP2Z
      IF(LT.EQ.0)GO TO 13
      WRITE(LUNOUT,50)HEF(J),PXF(J),PYF(J),PZF(J)
      WRITE(LUNOUT,50)HEF(I),PXF(I),PYF(I),PZF(I)
   50 FORMAT(1H0,11HE,PX,PY,PZ=,4F8.4)
   13 CONTINUE
      GO TO 200
  100 CONTINUE
C*****ONLY ONE RESONANCE WILL BE CREATED IF IT=0
      IF(LT.EQ.0)GO TO 14
      WRITE(LUNOUT,70)
   70 FORMAT(1H0,'CUT OFF OF THE LEFT AND RIGHT JET BEFOR THE FIRST
     *DECAY STEP',/,'IF ALLOWED ONLY ONE RESONANCE COULD BE CREATED')
   14 CONTINUE
      IHAD=1
*
*  La = 2 means resonance creation
*
      LA=2
  200 CONTINUE
      RETURN
      END
