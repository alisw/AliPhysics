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
*$ CREATE ANKEKA.FOR
*COPY ANKEKA
*                                                                      *
*=== ankeka ===========================================================*
*                                                                      *
      FUNCTION ANKEKA(IT,TO,AMSS,SQAMSS)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C     LAST CHANGE 05. DECEMBER 83    PERTTI AARNIO
C                                    HELSINKI UNIVERSITY OF
C                                    TECHNOLOGY, FINLAND
C
C     TO BE CALLED FROM THE HIGH ENERGY PRODUCTION
C
C     TO GIVE THE PROTON AND NEUTRON MULTIPLICITIES OF THE HIGH
C     ENERGY PART OF THE INTRANUCLEAR CASCADE  (N-2P,N-2N).
C     NOTE THAT IT IS SUPPOSED IN EKEKA AND RAKEKA THAT
C     N-1P=9*N-2P, N-1N=9*N-2N.
C
C     NOTE!!!!!!! NON-STANDARD PARTICLE NUMBERING
C
C     INPUT VARIABLES:
C        IT     = TYPE OF THE PARTICLE
C                 1=INTRANUCLEAR PROTON
C                 2=INTRANUCLEAR NEUTRON
C        TO     = ENERGY OF THE COLLIDING PARTICLE IN GEV
C        AMSS   = ATOMIC WEIGHT OF THE MEDIUM
C        SQAMSS = SQRT(AMSS)
C
C     SEE RANFT/ROUTTI PARTICLE ACC. VOL 4 P 105
C********************************************************************
C
      DIMENSION A(2),B(2),C(2),D(2)
      SAVE A,B,C,D
C
C
C     THESE ARE EFFECTIVE FLUKA81 MULTIPLICITIES
C
      DATA A/1.D0,1.3D0/
      DATA B/0.06D0,0.09D0/
      DATA C/1.104D0,1.90D0/
C     DATA D/0.14D0,0.19D0/
C  This ones for continuity !!
      DATA D/0.245333D0,0.33333D0/
      IF(TO.GT.0.1D0)GO TO 1
      ANKEKA=SQAMSS*B(IT)*0.1D0
      RETURN
 1    CONTINUE
      IF(TO.GE.10.D0)GO TO 2
      ANKEKA=0.1D0*SQAMSS*(0.5D0+A(IT)*(1.D0+LOG10(TO))**2)*D(IT)
      RETURN
 2    CONTINUE
      ANKEKA=SQAMSS*0.1D0*C(IT)
      RETURN
      END
