*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:48  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if !defined(CERNLIB_QMILSTD)
      FUNCTION ISHFTC(I,K,L)
C
C CERN PROGLIB# M441    ISHFTC          .VERSION KERNFOR  4.16  870601
C ORIG. 01/10/85  T.LINDELOF
C
C      FORTRAN VERSION OF MIL-STD-1753 ISHFTC FUNCTION
C      K IS THE SHIFT COUNT - NEGATIVE TO THE RIGHT
C                             POSITIVE TO THE LEFT
C      L IS THE NUMBER OF THE RIGHTMOST BITS TO BE SHIFTED, ALL
C         OTHER BITS ARE LEFT UNCHANGED
C
#include "kerngen/wordsize.inc"
#include "kerngen/iall1s.inc"

#include "kerngen/q_andor.inc"
* Ignoring t=pass
C
      IF (K.GT.0) THEN
         ISHFTC = IOR(IAND(ISHFT(IALL11, L), I),
     +    IOR(ISHFT(ISHFT(I, NBITPW-L), K-NBITPW),
     +        ISHFT(ISHFT(I, NBITPW-L+K), L-NBITPW)))
      ELSEIF (K.LT.0) THEN
         ISHFTC = IOR(IAND(ISHFT(IALL11, L), I),
     +        IOR(ISHFT(ISHFT(I, NBITPW-L), L+K-NBITPW),
     +        ISHFT(ISHFT(I, NBITPW+K), L-NBITPW) ))
      ELSE
         ISHFTC = I
      ENDIF
      END
#endif
