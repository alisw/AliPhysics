*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      SUBROUTINE GPFIS
C
C *** GENERATION OF PHOTO-FISSION AND PHOTO-ABSORBTION MECHANISMS ***
C *** HMF 25-AUG-1989 RWTH AACHEN / NVE 11-MAY-1990 CERN GENEVA ***
C
C CALLED BY : GTGAMA
C ORIGIN    : H.FESEFELDT (ROUTINE IPFISS)
C
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcking.inc"
#include "geant321/gckine.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gsecti.inc"
C
      DIMENSION RNDM(3),SPNEUT(10)
      LOGICAL CALFL
      SAVE SPNEUT,CALFL
C
      DATA SPNEUT/10*0./
      DATA CALFL/.FALSE./
C
      KCASE = NAMEC(23)
C
      IF(IPFIS.NE.1) THEN
        DESTEP = DESTEP + GETOT
        GOTO 90
      ENDIF
C
      ISTOP = 1
C
C     SELECT SUBPROCESS
C
      STEPAB = (1.-GEKRAT)* Q(JPFIS+NEK1+IEKBIN)
     +       +     GEKRAT * Q(JPFIS+NEK1+IEKBIN+1)
      CALL GRNDM(RNDM,1)
      IF(RNDM(1).LT.STEPPF/STEPAB) GOTO 10
C
C     PHOTOABSORBTION
C
      NGKINE    = 1
      TOFD(1)   = 0.0
      GKIN(5,1) = 13
      CALL GRNDM(RNDM,3)
      JPA       = LQ(JPART-13)
      GKIN(4,1) = Q(JPA+7) - 0.002*LOG(RNDM(1))
      COST      = -1.+2.*RNDM(2)
      SINT      = SQRT(1.-COST*COST)
      PHI       = TWOPI*RNDM(3)
      PPHMF     = GKIN(4,1)**2-Q(JPA+7)**2
      IF(PPHMF.LT.0.) PPHMF=0.
      PPHMF     = SQRT(PPHMF)
      GKIN(1,1) = PPHMF*SINT*SIN(PHI)
      GKIN(2,1) = PPHMF*SINT*COS(PHI)
      GKIN(3,1) = PPHMF*COST
      GPOS(1,1) = VECT(1)
      GPOS(2,1) = VECT(2)
      GPOS(3,1) = VECT(3)
C
      GOTO 100
C
C     PHOTOFISSION
C
   10 IF(CALFL) GOTO 20
      CALFL=.TRUE.
      XX        = 1.-0.5
      XXX       = SQRT(2.29*XX)
      SPNEUT(1) = EXP(-XX/0.965)*(EXP(XXX)-EXP(-XXX))/2.
      DO 11 I=2,10
      XX        = I*1.-0.5
      XXX       = SQRT(2.29*XX)
   11 SPNEUT(I) = SPNEUT(I-1)+EXP(-XX/0.965)*(EXP(XXX)-EXP(-XXX))/2.
      DO 12 I=1,10
   12 SPNEUT(I) = SPNEUT(I)/SPNEUT(10)
C
C
   20 NGKINE    = 1
      GKIN(1,1) = VECT(4)*VECT(7)
      GKIN(2,1) = VECT(5)*VECT(7)
      GKIN(3,1) = VECT(6)*VECT(7)
      GKIN(4,1) = GETOT
      GKIN(5,1) = 1
      TOFD(1)   = 0.0
      GPOS(1,1) = VECT(1)
      GPOS(2,1) = VECT(2)
      GPOS(3,1) = VECT(3)
C
C     NUMBER OF NEUTRONS AND PHOTONS
C
      EKHMF     = GEKIN*1000.
      IF ( EKHMF.LT.1.) EKHMF=1.
      AVERN     = 2.569+0.559*LOG(EKHMF)
      AVERG     = 0.500+0.600*LOG(EKHMF)
      CALL NORMAL(RAN)
      NN        = IFIX(AVERN+RAN*1.23+0.5)
      CALL NORMAL(RAN)
      NG        = IFIX(AVERG+RAN*3.00+0.5)
      IF(NN.LT.1) NN=1
      IF(NG.LT.1) NG=1
C
C     DISTRIBUTE KINETIC ENERGY
C
      JPA = LQ(JPART-13)
C
      DO 25 I=1,NN
C --- Protect against stack overflow ---
      IF (NGKINE .GE. MXGKIN) GO TO 31
         NGKINE = NGKINE+1
         CALL GRNDM(RNDM,1)
         DO 21 J=1,10
            IF(RNDM(1).LT.SPNEUT(J)) GOTO 22
   21    CONTINUE
         J=10
   22    CALL GRNDM(RNDM,1)
         GKIN(4,NGKINE) = (J-1)*1. + RNDM(1) + Q(JPA+7)*1000.
         GKIN(5,NGKINE) = 13
         TOFD(NGKINE)   = 0.0
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
   25 CONTINUE
C
      DO 30 I=1,NG
C --- Protect against stack overflow ---
      IF (NGKINE .GE. MXGKIN) GO TO 31
         NGKINE=NGKINE+1
         CALL GRNDM(RNDM,1)
         GKIN(4,NGKINE) = -0.87*LOG(RNDM(1))
         GKIN(5,NGKINE) = 1
         CALL GRNDM(RNDM,1)
         TOFD(NGKINE)   = -25.E-9*LOG(RNDM(1))
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
   30 CONTINUE
C
C --- GO BACK TO GEV UNITS ---
 31   CONTINUE
      DO 35 I=2,NGKINE
         GKIN(4,I)=GKIN(4,I)*1E-3
   35 CONTINUE
C
C     DISTRIBUTE DIRECTIONS ISOTROPICALLY IN LAB- SYSTEM
C
      DO 40 I=2,NGKINE
         CALL GRNDM(RNDM,1)
         COST      = -1.+2.*RNDM(1)
         SINT      = SQRT(1.-COST*COST)
         CALL GRNDM(RNDM,1)
         PHI       = TWOPI*RNDM(1)
         IF (GKIN(5,I).LT.1.5) THEN
            PPHMF  = GKIN(4,I)
         ELSE
            PPHMF  = GKIN(4,I)**2 - Q(JPA+7)**2
            IF(PPHMF.LT.0.) PPHMF=0.
            PPHMF  = SQRT(PPHMF)
         ENDIF
         GKIN(1,I) = PPHMF*SINT*SIN(PHI)
         GKIN(2,I) = PPHMF*SINT*COS(PHI)
         GKIN(3,I) = PPHMF*COST
   40 CONTINUE
C
      GOTO 100
C
   90 ISTOP = 2
      IPART=1
      JPA = LQ(JPART-IPART)
      DO 95 J=1,5
         NAPART(J) = IQ(JPA+J)
   95 CONTINUE
      ITRTYP = Q(JPA+6)
      AMASS  = Q(JPA+7)
      CHARGE = Q(JPA+8)
      TLIFE  = Q(JPA+9)
      VECT(7)= 0.0
      GETOT  = 0.0
      GEKIN  = 0.0
C
 100  CONTINUE
      END
