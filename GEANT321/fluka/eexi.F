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
*$ CREATE EEXI.FOR
*COPY EEXI
*                                                                      *
*=== eexi =============================================================*
*                                                                      *
      FUNCTION EEXI(IJ,E,A)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C     VERSION BY                     J. RANFT
C                                    LEIPZIG
C
C     NUCLEAR EXCITATION ENERGY INCLUDING ALL LOW ENERGY SECONDARIES
C     IN COLLISIONS OF PARTICLE IJ OF ENERGY E GEV ON NUCLEI A
C
C     NOTE!!!!!!! NON-STANDARD PARTICLE NUMBERING
C
C********************************************************************
C
      IF (E.LE.0.125D0) GO TO 20
      B=SQRT(A)/9.D0-.2D0
      IF (B.GT.A*0.01D0) B=A*0.01D0
      IF (E.GE.3.D0) GO TO 10
      IF (B.LT.0.125D0) B=0.125D0
      EEXI=0.125D0+(E-0.125D0)*(B-0.125D0)/2.875D0
      RETURN
   10 CONTINUE
      EEXI=B
      RETURN
   20 CONTINUE
      EEXI=E
      RETURN
      END
