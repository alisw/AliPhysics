*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMDOS)
#include "dosgs/ibits.F"
#elif defined(CERNLIB_QMLNX)
#include "lnxgs/ibits.F"
#elif !defined(CERNLIB_QMILSTD)
      FUNCTION IBITS(I,N,L)
C
C CERN PROGLIB# M441    IBITS           .VERSION KERNFOR  4.16  870601
C ORIG. 01/03/85  M.METCALF
C
C      FORTRAN VERSION OF MIL-STD-1753 FUNCTION IBITS
C      N IS THE RIGHTMOST BIT OF THE REQUIRED FIELD
C      L IS THE LENGTH OF THE REQUIRED FIELD
C
#include "kerngen/wordsize.inc"
C
      IBITS = ISHFT(ISHFT(I, NBITPW-L-N), L-NBITPW)
      END
#endif
