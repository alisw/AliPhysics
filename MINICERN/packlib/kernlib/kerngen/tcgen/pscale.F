*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      FUNCTION PSCALE (NSC,NMAX,A,NST)
C
C CERN PROGLIB# M215    PSCALE          .VERSION KERNFOR  1.0   681105
C ORIG. 05/11/68
C
C--                FIND THE LARGEST POWER PSCALE=10.**NSC
C-                 SUCH THAT A*PSCALE HAS AT MOST NST DIGITS ON THE
C-                                    LEFT OF THE DECIMAL POINT
C-                 NSC LIMITED TO NSC.LE.NMAX
C
C
      NSC=NMAX
      IF (A.EQ.0.) GO TO 9
      N=100.+ LOG10 (ABS (A))+.0001
      N=NST+99-N
      IF (N.LT.NSC) NSC=N
    9 PSCALE=10.**NSC
      RETURN
      END
