*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_QMILSTD)
      SUBROUTINE MVBITS(I,N,L,J,M)
C
C CERN PROGLIB# M441    MVBITS          .VERSION KERNFOR  4.16  870601
C ORIG. 01/10/85  T.LINDELOF
C
C      FORTRAN VERSION OF MIL-STD-1753 SUBROUTINE MVBITS
C      N IS THE POSITION IN I FROM WHICH BITS ARE EXTRACTED
C      L IS THE NUMBER OF BITS TO BE EXTRACTED
C      M IS THE POSITION IN J INTO WHICH BITS ARE TO BE PLACED
C
#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"

#include "kerngen/q_andor.inc"
* Ignoring t=pass
C
      J = IOR(IAND(ISHFT(IALL11, L+M), J),
     +       IOR(IAND(NOT(ISHFT(IALL11, M)), J),
     +       ISHFT(ISHFT(ISHFT(I, NBITPW-L-N), L-NBITPW), M) ))
      END
#endif
