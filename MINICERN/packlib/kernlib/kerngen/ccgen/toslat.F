*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:28  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
      SUBROUTINE TOSLAT (IVEC,N)
C
C CERN PROGLIB#         TOSLAT          .VERSION KERNFOR  4.33  920511
C ORIG. 22/02/91, JZ
C Fortran interface to copy a vector from a C routine into /SLATE/
C On Apollo it is at the moment not possible to access a
C Fortran COMMON from C
C
      COMMON /SLATE/ISLAT(40)

      CALL UCOPY (IVEC, ISLAT, N)
      RETURN
      END
