*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:22  mclareni
* Kernlib
*
*
      FUNCTION RNDM (ISEED)
C
C CERN PROGLIB# V104    RNDM            .VERSION KERNIRT  1.06  930811
C ORIG. 22-MAR-88 from IBM version of G.Marsaglia et al., Montreal
C
C     !!!   Version valid only for 32 bit machines   !!!
C     !!!   with IEEE floating point representation  !!!
C
C-    Uniform Random Number Generator,
C-    giving the same sequence as the IBM and VAX version

      EQUIVALENCE (AMAN,MANT)
      PARAMETER    (MSK1 =  25165824)
      PARAMETER    (MSK2 = 855638016)
C     PARAMETER    (MSK1 = '0C000000'X)
C     PARAMETER    (MSK2 = '33000000'X)
      SAVE  MCGN
      DATA  MCGN  /12345/

#if !defined(CERNLIB_QISASTD)
#include "kerngen/q_andor.inc"
#include "kerngen/q_shift.inc"
* Ignoring t=pass
#endif

      MCGN = MCGN * 69069
#if defined(CERNLIB_QISASTD)
      MANT = ISHFT (MCGN,-8)
#endif
#if !defined(CERNLIB_QISASTD)
      MANT = ISHFTR(MCGN, 8)
#endif
      IF (MANT.EQ.0)         GO TO 14
      AMAN = MANT
C-    AMAN in the range 1 to 2**24-1
      MANT = MANT - MSK1
C-    multiply by 2.**(-24)
      RNDM = AMAN
      RETURN

C--   for zero set RNDM = 2.**(-25)
   14 RNDM =  0.29802322387695312E-07

      RETURN

C--       Integer random number
      ENTRY IRNDM (ISEED)
      MCGN  = MCGN * 69069
#if defined(CERNLIB_QISASTD)
      IRNDM = ISHFT (MCGN,-1)
#endif
#if !defined(CERNLIB_QISASTD)
      IRNDM = ISHFTR(MCGN, 1)
#endif
      RETURN

C--       Set the seed
      ENTRY RDMIN (ISEED)
      MCGN = ISEED
      RETURN

C--       Get the seed
      ENTRY RDMOUT (ISEED)
      ISEED = MCGN
      END
