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
*$ CREATE EKEKA.FOR
*COPY EKEKA
*                                                                      *
*=== ekeka ============================================================*
*                                                                      *
      FUNCTION EKEKA(IX,TO,AMSS,SQAMSS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C     LAST CHANGE 05. DECEMBER BY    PERTTI AARNIO
C                                    HELSINKI UNIVERSITY OF
C                                    TECHNOLOGY, FINLAND
C
C     TO BE CALLED FROM THE HIGH ENERGY PRODUCTION
C
C     THIS IS A SUBROUTINE OF FLUKA TO CALCULATE THE ENERGY AVAILABLE
C     FOR THE REACTION "IX".
C
C     NOTE!!!!!!! REACTION NUMBERING IS NOT
C     COMPATIBLE WITH PARTICLE NUMBERING
C
C     INPUT VARIABLES:
C        IX     = TYPE OF THE REACTION
C                 1=NUCLEAR EXCITATION
C                 2=INTRANUCLEAR PROTON
C                 3=INTRANUCLEAR NEUTRON
C                 4=1+2
C                 5=1+2+3
C        TO     = KINETIC ENERGY OF THE COLLIDING PARTICLE IN GEV
C        AMSS   = ATOMIC WEIGHT OF THE MEDIUM
C        SQAMSS = SQRT(AMSS)
C
C     SEE RANFT/ROUTTI PARTICLE ACC VOL 4 P 106
C
C     NOTE THAT IN INTRANUCLEAR PART AVERAGE TOTAL ENERGY IS
C     OBTAINED BY MULTPLYING THE AVERAGE ENERGY OF THE HIGH ENERGY
C     PARTICLES BY THE MULTIPLICITY OF THE HIGH ENERGY PARTICLES
C     I.1. E-TOT,AV=2.5*N2*E-AV(ALFA-2).  THE FACTOR 2.5 IS
C     NEEDED TO TAKE INTO ACCOUNT THE LOW ENERGY PART ALSO.
C     2.5 IS BASED ON THE ASSUMPTION THAT N1/N2=9 AND THAT
C     (E-AV(ALFA-2))/(E-AV(ALFA-1))=6.
C********************************************************************
C
      GO TO (1,2,3,1,1),IX
C
C
    1 CONTINUE
      IF (TO.GT.01D0)GO TO 12
      AA=0.001D0*SQAMSS
      GO TO 19
   12 CONTINUE
      APAR=0.035D0
      BPAR=3.D0
      CPAR=0.1D0
      AA=CPAR*SQAMSS*(0.01D0+APAR*(BPAR+LOG10(TO))**2)
C
   19 CONTINUE
      IF (IX.GT.3) GO TO 2
      EKEKA=AA
      RETURN
C
C
    2 CONTINUE
      AN=ANKEKA(1,TO,AMSS,SQAMSS)
      A=AKEKA(1,TO,AMSS)
      EXTOA=0.D0
      IF(TO.LT.5.D0*A) EXTOA=EXP(-TO/A)
      TPKAV=A*(1.D0-(TO/A+1.)*EXTOA)/(1.D0-EXTOA)
      BB=2.5D0*TPKAV*AN
C
      IF (IX.EQ.4) GO TO 4
      IF (IX.EQ.5) GO TO 3
      EKEKA=BB
      RETURN
C
C
    3 CONTINUE
      AN=ANKEKA(2,TO,AMSS,SQAMSS)
      A=AKEKA(2,TO,AMSS)
      EXTOA=0.D0
      IF (TO.LT.5.D0*A) EXTOA=EXP(-TO/A)
      TNKAV=A*(1.D0-(TO/A+1.D0)*EXTOA)/(1.D0-EXTOA)
      CC=2.5D0*TNKAV*AN
C
      IF (IX.EQ.5) GO TO 5
      EKEKA=CC
      RETURN
C
C
    4 CONTINUE
      EKEKA=AA+BB
      RETURN
C
C
    5 CONTINUE
      EKEKA=AA+BB+CC
      RETURN
      END
