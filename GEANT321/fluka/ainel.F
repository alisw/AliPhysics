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
*$ CREATE AINEL.FOR
*COPY AINEL
*                                                                      *
*=== ainel ============================================================*
*                                                                      *
      FUNCTION AINEL(IJ,ILO,E,A,SQA)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C     LAST CHANGE 16. JULY 81 BY     PERTTI AARNIO
C                                    HELSINKI UNIVERSITY OF
C                                    TECHNOLOGY, FINLAND
C
C
C     FUNCTION OF FLUKA TO GIVE INELASTICITIES FOR PRODUCTION OF
C     PARTICLES ILO BY PARTICLE IJ
C
C     NOTE!!!!!!! NON-STANDARD PARTICLE NUMBERING
C
C     INPUT VARIABLES:
C        IJ    = TYPE OF THE PRIMARY
C                1=PROTON
C                2=NEUTRON
C                3=CHARGED PION
C
C        ILO   = TYPE OF THE SECONDARY
C                1=PROTON
C                2=NEUTRON
C                3=CHARGED PION
C                4=NEUTRAL PION
C                5=NUCLEAR EXCITATION EXCLUDING LOW ENERGY
C                  SECONDARIES FROM INTRANUCLEAR CASCADE
C                6=INTRANUCLEAR PROTON
C                7=INTRANUCLEAR NEUTRON
C
C        E     = KINETIC ENERGY OF THE PRIMARY IN GEV
C        A     = ATOMIC WEIGHT OF THE MEDIUM
C        SQA   = SQRT(A)
C
C     OTHER VARIABLES:
C        AK    = LOOK-UP TABLE FOR INELASTCITIES
C
C********************************************************************
C
      DIMENSION AK(4,3)
      SAVE AK
      DATA AK/
     10.35D0,0.3D0,0.25D0,0.1D0,
     10.3D0,0.35D0,0.25D0,0.1D0,
     10.15D0,0.15D0,0.5D0,0.2D0/
C
C
      GO TO (1234,1234,1234,1234,5,6,7),ILO
C
C
 1234 CONTINUE
      IF(E.LT.0.125D0) GO TO 10
      AINEL=AK(ILO,IJ)*(1.D0-EEXI(IJ,E,A)/E)
      RETURN
C
C
 5    CONTINUE
      IF(E.LT.0.125D0) GO TO 11
      ANEL=EEXI(IJ,E,A)/E
      AKEK=EKEKA(2,E,A,SQA)/E+EKEKA(3,E,A,SQA)/E
      AINEL=ANEL-AKEK
      IF (AKEK.GE.ANEL) AINEL=EKEKA(1,E,A,SQA)/E
      RETURN
C
C
    6 CONTINUE
      IF(E.LT.0.125D0) GO TO 10
      AINEL=EKEKA(2,E,A,SQA)/E
      RETURN
C
C
    7 CONTINUE
      IF(E.LT.0.125D0) GO TO 10
      AINEL=EKEKA(3,E,A,SQA)/E
      RETURN
C
C
   10 CONTINUE
      AINEL=0.D0
      RETURN
C
C
 11   CONTINUE
      AINEL=1.D0
      RETURN
      END
