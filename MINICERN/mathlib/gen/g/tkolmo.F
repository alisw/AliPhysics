*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:42  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE TKOLMO(A,NA,B,NB,PROB)
*==========>
*        Tests whether two one-dimensional sets of points
*        are compatible with coming from the same parent
*        distribution, using the Kolmogorov test
*   A is the vector of NA points in set one.
*   B is the vector of NB points in set two.
*   PROB is the probability of compatibility returned to
*      the calling program.
*
*   F. James, April, 1987
*==========>
*
      DIMENSION A(NA),B(NB)
*                 In case of an error, PROB is returned as -1.
      PROB = -1.0
*                 Require at least two points in each set
      IF (NA.LE.2 .OR. NB.LE.2)  GO TO 99
*                 Constants needed
      RNA = REAL(NA)
      RNB = REAL(NB)
      SA = 1.0/RNA
      SB = 1.0/RNB
*                 Starting values for main loop
      IF (A(1) .LT. B(1))  THEN
         RDIFF = -SA
         IA = 2
         IB = 1
      ELSE
         RDIFF = SB
         IB = 2
         IA = 1
      ENDIF
      RDMAX = ABS(RDIFF)
*
*             Main loop over point sets to find max distance
*             RDIFF is the running difference, and RDMAX the max.
      DO 10 I= 1, NA+NB
      IF (A(IA) .LT. B(IB)) THEN
         RDIFF = RDIFF - SA
         IA = IA + 1
         IF (IA .GT. NA)  GO TO 20
      ELSE
         RDIFF = RDIFF + SB
         IB = IB + 1
         IF (IB .GT. NB)  GO TO 20
      ENDIF
      RDMAX = MAX(RDMAX,ABS(RDIFF))
   10 CONTINUE
*            Should never terminate this loop!
      GO TO 99
   20 RDMAX = MAX(RDMAX,ABS(RDIFF))
      Z = RDMAX * SQRT(RNA*RNB/(RNA+RNB))
      PROB = PROBKL(Z)
*
   99 CONTINUE
      END
