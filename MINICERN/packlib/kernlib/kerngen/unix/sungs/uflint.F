*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:08  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE UFLINT (VECT,NWP,MODEP)
C
C CERN PROGLIB# M251    UFLINT          .VERSION KERNSUN  1.00  880922
C ORIG. 22/09/88, JZ
C
C-    CONVERT VECT(1->NWP) IN SITU TO FLOATING OR INTEGERS ACCORDING
C-    TO MODEP
C
      DIMENSION    VECT(99), NWP(9), MODEP(9)
      EQUIVALENCE  (A,IA), (B,IB)
      LOGICAL      VARY
C
C
      NW   = -NWP(1)
      MODE = MODEP(1)
      N    = 0
      VARY = .TRUE.
      IF (NW.GE.0)           GO TO 31
      NW   = -NW
      VARY = .FALSE.
      JUSE = AND (MODE,1)
      GO TO 31
C
   30 VECT(N) = B
   31 IF (N.EQ.NW)           RETURN
      N = N + 1
      IF (.NOT.VARY)         GO TO 34
      JUSE = AND (MODE,1)
      MODE = RSHIFT (MODE,1)
   34 A    = VECT(N)
      JCH  = RSHIFT (IA,23)
      IF (JCH.EQ.0)          GO TO 36
      IF (JCH.EQ.-1)         GO TO 36
C
C--                CONVERT TO INTEGER
C
      IF (JUSE.EQ.0)         GO TO 31
      IB = A
      GO TO 30
C
C--                CONVERT TO FLOATING
   36 IF (JUSE.NE.0)         GO TO 31
      B = IA
      GO TO 30
      END
