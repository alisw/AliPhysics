*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
*
      SUBROUTINE GVPSQR (CBA, SOL, NSOL)
************************************************************************
*                                                                      *
*     SUBROUTINE GVPSQR (CBA,SOL,NSOL)             870924  VP          *
*                                                                      *
*       SOLVE  QUADRATIC  EQUATION                                     *
*                                                                      *
*   ARGUMENTS:                                                         *
*       CBA     Array of coeff's A0 + A1*x + A2*x**2                   *
*       SOL     Solutions                                              *
*       NSOL    Number of solutions :                                  *
*               if zero - SOL(1)= extremum                             *
*               if -ve  - No solution at all                           *
*                                                                      *
************************************************************************
 
      PARAMETER      (ZERO2=1.E-12)
      REAL           CBA(3), SOL(2), SWAP
*-----------------------------------------------------------------------
 
      A = CBA(3)
      B = CBA(2)*.5
      C = CBA(1)
      IF (B.LT.0.)              THEN
        A = -A
        B = -B
        C = -C
      ENDIF
 
      AMX = MAX (ABS(A),ABS(B),ABS(C))
      IF (AMX.LE.0.)            GO TO 100
      A   = A/AMX
      B   = B/AMX
      C   = C/AMX
 
      DIS = B**2 - A*C
      NSOL = 1
      IF (ABS(DIS) .LE. ZERO2)  DIS = 0.
      IF (DIS.LT.0.)            THEN
        NSOL = 0
        DIS  = 0.
      ENDIF
 
      IF (DIS.GT.0.)            DIS = SQRT(DIS)
      BDIS = B + DIS
      IF (ABS(C).GT.1.E+10*BDIS)        GO TO 100
      SOL(1) = 0.
      IF (ABS(BDIS).LE.0.)      GO TO 999
      SOL(1) = (-C/BDIS)
      IF (DIS.LE.0.)            GO TO 999
      IF (BDIS.GE.1.E+10*ABS(A))        GO TO 999
      NSOL   = 2
      SOL(2) = (-BDIS/A)
      IF (SOL(1).GT.SOL(2))     THEN
        SWAP = SOL(1)
        SOL(1) = SOL(2)
        SOL(2) = SWAP
      ENDIF
      GO TO 999
 
*
**      No solutions
  100 NSOL = -1
 
  999 RETURN
      END
