*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:06  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*=== flavor ===========================================================*
*
      SUBROUTINE FKFLAV(IT,LT,LL,E0,IV,RE,KFR1,KFR2,ISU,BETE,KFA1,KFA2,
     &                  KFA3,KFA4,IOPT)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     New version from A. Ferrari ( 22 august 1990 ): it is almost the *
*     the same as before, but with a few corrections important for isu *
*     2 and 3 and also for isu=4 if it switchs to 100 continue         *
*                                                                      *
*----------------------------------------------------------------------*
#include "geant321/part.inc"
*
*   The following are the masses of the quarks: the d quark mass is
*   assumed to be the same as the u one. They are quite different from
*   the last values from the particle data group, but any change
*   can imply a change also in the BET parameter in common INPDAT
      PARAMETER ( UQUARM = 0.3D+00 )
      PARAMETER ( SQUARM = 0.5D+00 )
      PARAMETER ( CQUARM = 2.1D+00 )
      PARAMETER ( BQUARM = 5.0D+00 )
*
      DIMENSION RE(*),KFR1(*),KFR2(*),IV(*)
      REAL RNDM(2)
*
C     CHOICE OF THE QUARK FLAVOUR
      IF (LT.EQ.1)WRITE(LUNOUT,288)IT,LT,LL,E0,ISU,BETE,KFA1,KFA2
  288 FORMAT(3I5,E12.4,I5,E12.4,2I5
     *,' FLAVOR IT,LT,LL,E0,ISU,BETE,KFA1,KFA2')
      I=IT
      J=IT-1
      IVA=1
      IVX=IV(I)
      IF (I .LE. 1) THEN
         IF (IOPT.EQ.2) THEN
            KX1 = KFA1
            KX2 = KFA2
         ELSE IF (IOPT.EQ.3.AND.LL.EQ.1) THEN
            KX1=KFA2
            KX2=0
         ELSE IF (IOPT.EQ.4 .AND. KFA1.LE.6 .AND. LL.EQ.1) THEN
            KX1=KFA2
            KX2=KFA3
         ELSE IF (IOPT.EQ.4 .AND. KFA1.GT.6 .AND. LL.EQ.0) THEN
            KX1=KFA2
            KX2=KFA3
         ELSE IF (IOPT.EQ.5 .AND. LL.EQ.0) THEN
            KX1=KFA3
            KX2=KFA4
         ELSE IF (IOPT.EQ.5 .AND. LL.EQ.1) THEN
            KX1=KFA1
            KX2=KFA2
         ELSE
            KX1=KFA1
            KX2=0
         END IF
         RX = E0
      ELSE
         KX1=KFR1(J)
         KX2=KFR2(J)
         RX =RE(J)
      END IF
      IF (KX1.GT.0.AND.KX2.GT.0) THEN
         BET=10.D+00
      ELSE
         BET=BETE
      END IF
      CALL GRNDM(RNDM,2)
      Z1=RNDM(1)
      Z2=RNDM(2)
      IF(ISU.EQ.4) GO TO 300
      IF(ISU.EQ.3) GO TO 200
      IF(ISU.EQ.2) GO TO 100
C     U FLAVOUR
         KF1=1
         KF2=1
      GO TO 20
C     U/D FLAVOURS
  100 CONTINUE
         IF (KX1.EQ.1.OR.KX1.EQ.7) THEN
            IIAA=1
         ELSE IF (KX1.EQ.2.OR.KX1.EQ.8) THEN
            IIAA=2
         ELSE
            IIAA=0
         END IF
  110 CONTINUE
         IF (IIAA .EQ. 1) THEN
            PD=0.6666666666666667D+00
            PU=0.3333333333333333D+00
         ELSE IF (IIAA .EQ. 2) THEN
            PU=0.6666666666666667D+00
            PD=0.3333333333333333D+00
         ELSE
            PU=0.5D+00
            PD=0.5D+00
         END IF
         PS=0.D+00
         PC=0.D+00
         IF (Z1 .LE. PD) THEN
            KF1=2
         ELSE
            KF1=1
         END IF
         IF (Z2 .LE. PD) THEN
            KF2=2
         ELSE
            KF2=1
         END IF
      GO TO 20
C     U/D/S FLAVOURS
  200 CONTINUE
         IF (KX1.EQ.1.OR.KX1.EQ.7) THEN
            IIAA=1
         ELSE IF (KX1.EQ.2.OR.KX1.EQ.8) THEN
            IIAA=2
         ELSE
            IIAA=0
         END IF
  210 CONTINUE
         IF (RX .LE. 1.019D+00 ) GO TO 110
         X1=RX
         X2=UQUARM
         PU=BETA(X1,X2,BET)
         X2=SQUARM
         PS=BETA(X1,X2,BET)
         PTOT=2.D+00*PU+PS
         PU1=PU/PTOT
         PS =PS/PTOT
         PC=0.D+00
         IF (IIAA .EQ. 1) THEN
            PU=0.6666666666666667D+00*PU1
            PD=2.D+00*PU1-PU
         ELSE IF (IIAA .EQ. 2) THEN
            PD=0.6666666666666667D+00*PU1
            PU=2.D+00*PU1-PD
         ELSE
            PU=PU1
            PD=PU
         END IF
         IF (Z1 .LE. PU) THEN
            KF1 = 1
         ELSE IF ( Z1 .LE. PU + PD ) THEN
            KF1 = 2
         ELSE
            KF1 = 3
         END IF
         IF (Z2 .LE. PU) THEN
            KF2 = 1
         ELSE IF ( Z2 .LE. PU + PD ) THEN
            KF2 = 2
         ELSE
            KF2 = 3
         END IF
      GO TO 20
C     U/D/S/C FLAVOUR
  300 CONTINUE
         GO TO (11,12,13,14,14,11,12,13,14,14),IVX
   11    CONTINUE
            IF (KX1.EQ.4.OR.KX1.EQ.10) THEN
               GM=AM(129)
            ELSE
               GM=AM(127)
            END IF
            IF (KX1.EQ.1.OR.KX1.EQ.7) THEN
               IIAA=1
            ELSE IF (KX1.EQ.2.OR.KX1.EQ.8) THEN
               IIAA=2
            ELSE
               IIAA=0
            END IF
         GO TO 15
   12    CONTINUE
            IF (KX1.EQ.4.OR.KX1.EQ.10) THEN
               GM=AM(170)
            ELSE
               GM=AM(127)
            END IF
            IIAA=0
         GO TO 15
   13    CONTINUE
*  |  |  +-------------------------------------------------------------*
*  |  |  | The following if replaces the cards:
*  |  |  |      GM=3.85D0
*  |  |  |      IF(KX1.EQ.4.AND.KX2.EQ.4) GM=4.89D0
*  |  |  |      IF(KX1.NE.4.AND.KX2.NE.4) GM=2.770D0
*  |  |  |      IF(KX1.EQ.10.AND.KX2.EQ.10) GM=4.89D0
*  |  |  |      IF(KX1.NE.10.AND.KX2.NE.10) GM=2.770D0
*  |  |  | It is completely equivalent except for the combination
*  |  |  | 4-4 which now gives GM = 4.89, while in the original
*  |  |  | coding gave GM = 2.77, because of the last condition
*  |  |  | always overrides the first one (it seems to be a mistake)
            IF (KX1 .EQ. 4 .OR. KX1 .EQ. 10 ) THEN
               IF ( KX2 .EQ. KX1 ) THEN
                  GM = AM(170)
               ELSE
                  GM = AM(169)
               END IF
            ELSE IF (KX2 .EQ. 4 .OR. KX2 .EQ. 10 ) THEN
               GM = AM(169)
            ELSE
               GM = AM(166)
            END IF
            IIAA=0
         GO TO 15
   14    CONTINUE
*  |  |  +-------------------------------------------------------------*
*  |  |  | The following if replaces the cards:
*  |  |  |      GM=3.684D0
*  |  |  |      IF(KX1.NE.4.AND.KX2.NE.4) GM=2.140D0
*  |  |  |      IF(KX1.NE.10.AND.KX2.NE.10) GM=2.140D0
*  |  |  | It is equivalent: only for the combinations
*  |  |  | 4-10 and 10-4 we get GM = 3.684
*  |  |  | It is not clear if it is correct since 4-x,x-4 (x.ne.10),
*  |  |  | 10-x,x-10 (x.ne.4) give GM = 2.14
            IF ((KX1.EQ.4.AND.KX2.EQ.10).OR.(KX1.EQ.10.AND.KX2.EQ.4))
     &         THEN
               GM = AM(129)
            ELSE
               GM = AM(127)
            END IF
            IF (IVX.EQ.4.OR.IVX.EQ.9) THEN
               KAXI=KX1
            ELSE IF (IVX.EQ.5.OR.IVX.EQ.10) THEN
               KAXI=KX2
            ELSE
*  |  |  |   Kaxi = 0 added for completeness, maybe it is useless
               KAXI=0
            END IF
            IF (KAXI.EQ.1.OR.KAXI.EQ.7) THEN
               IIAA=1
            ELSE IF (KAXI.EQ.2.OR.KAXI.EQ.8) THEN
               IIAA=2
            ELSE
               IIAA=0
            END IF
         GO TO 15
   15    CONTINUE
         IF (RX .LE. GM) GO TO 200
         X1=RX
         X2=UQUARM
         PU=BETA(X1,X2,BET)
         X2=SQUARM
         PS=BETA(X1,X2,BET)
         X2=CQUARM
         PC=BETA(X1,X2,BET)
         PTOT=2.D+00*PU+PS+PC
         PU1=PU/PTOT
         PS=PS/PTOT
         PC=PC/PTOT
         IF (IIAA .EQ. 1) THEN
            PU=0.6666666666666667D+00*PU1
            PD=2.D+00*PU1-PU
         ELSE IF (IIAA .EQ. 2) THEN
            PD=0.6666666666666667D+00*PU1
            PU=2.D+00*PU1-PD
         ELSE
            PU=PU1
            PD=PU
         END IF
         IF (Z1 .LE. PU) THEN
            KF1 = 1
         ELSE IF ( Z1 .LE. PU + PD ) THEN
            KF1 = 2
         ELSE IF ( Z1 .LE. PU + PD + PS ) THEN
            KF1 = 3
         ELSE
            KF1 = 4
         END IF
         IF (Z2 .LE. PU) THEN
            KF2 = 1
         ELSE IF ( Z2 .LE. PU + PD ) THEN
            KF2 = 2
         ELSE IF ( Z2 .LE. PU + PD + PS ) THEN
            KF2 = 3
         ELSE
            KF2 = 4
         END IF
      GO TO 20
   20 CONTINUE
C*****CHOICE OF THE QUARKFLAVOURS IN DEPENDENCE OF THE VERTEX IV
      IVX=IV(I)
      GO TO (1,2,3,4,5,1,2,3,4,5),IVX
    1 CONTINUE
         IF (LL.EQ.1) THEN
            KFR1(I)=KF1+6
         ELSE
            KFR1(I)=KF1
         END IF
         KFR2(I)=0
      GO TO 30
    2 CONTINUE
         IF (LL.EQ.1) THEN
            KFR1(I)=KF1
            KFR2(I)=KF2
         ELSE
            KFR1(I)=KF1+6
            KFR2(I)=KF2+6
         END IF
      GO TO 30
    3 CONTINUE
         KFR2(I)=0
         IF (LL.EQ.1) THEN
            KFR1(I)=KF1+6
         ELSE
            KFR1(I)=KF1
         END IF
      GO TO 30
    4 CONTINUE
         IF (LL.EQ.1) THEN
            KFR1(I)=KF1
         ELSE
            KFR1(I)=KF1+6
         END IF
         KFR2(I)=KX2
      GO TO 30
    5 CONTINUE
         KFR1(I)=KX1
         IF (LL.EQ.1) THEN
            KFR2(I)=KF2
         ELSE
            KFR2(I)=KF2+6
         END IF
      GO TO 30
   30 CONTINUE
      IF(LT.EQ.0) GO TO 80
      WRITE(LUNOUT,6)PU,PD,PS,PC,KX1,KX2
    6 FORMAT(1H0,' FLAVOR PU,PD,PS,PC,KX1,KX2',4F8.4,2I5)
      IF(I.EQ.1) GO TO 40
      WRITE(LUNOUT,60)IV(I),LL,KFR1(J),KFR2(J),KFR1(I),KFR2(I)
      GO TO 50
   40 WRITE(LUNOUT,70)IV(I),LL,KFA1,KFA2,KFR1(I),KFR2(I)
   50 CONTINUE
   60 FORMAT(1H0,22HIV,LL,Q1A,Q2A,Q1N,Q2N=,6I3)
   70 FORMAT(1H0,'IV(I),LL,KFA1,KFA2,KFR1(I),KFR2(I)=',6I3 )
   80 CONTINUE
      RETURN
      END
