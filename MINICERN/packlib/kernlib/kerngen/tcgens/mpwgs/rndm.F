*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:24  mclareni
* Kernlib
*
*
      FUNCTION RNDM (ISEED)
C
C CERN PROGLIB# V104    RNDM            .VERSION KERNMPW  1.00  891208
C ORIG. 22-MAR-88 from IBM version of G.Marsaglia et al., Montreal
C
C-    Uniform Random Number Generator,
C-    giving the same sequence as the IBM and VAX version
 
      EQUIVALENCE (AMAN,MANT)
      PARAMETER    (MSK1 = 201326592)
      PARAMETER    (MSK2 = 855638016)
C     PARAMETER    (MSK1 = '0C000000'X)
C     PARAMETER    (MSK2 = '33000000'X)
      COMMON / V104RN / MCGN
      DATA MCGN  /12345/
 
 
      MCGN = MCGN * 69069
      MANT = ISHFT (MCGN,-8)
      IF (MANT.EQ.0)         GO TO 14
      AMAN = MANT
C-    AMAN in the range 1 to 2**24-1
      MANT = MANT - MSK1
C-    multiply by 2.**(-24)
      RNDM = AMAN
      RETURN
 
C--   for zero set RNDM = 2.**(-25)
   14 MANT = MSK2
      RNDM = AMAN
      RETURN
 
C--       Integer random number
      ENTRY IRNDM (ISEED)
      MCGN  = MCGN * 69069
      IRNDM = ISHFT (MCGN,-1)
      END
 
      SUBROUTINE RDMIN (ISEED)
C--       Set the seed
      COMMON / V104RN / MCGN
      MCGN = ISEED
      RETURN
 
C--       Get the seed
      ENTRY RDMOUT (ISEED)
      ISEED = MCGN
      END
