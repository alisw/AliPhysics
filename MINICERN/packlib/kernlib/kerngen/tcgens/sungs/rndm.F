*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:21  mclareni
* Kernlib
*
*
      FUNCTION RNDM (ISEED)
C
C CERN PROGLIB# V104    RNDM            .VERSION KERNSUN  1.06  920511
C ORIG. 13-JULY-88 JZ from IBM version of G.Marsaglia et al., Montreal
C
C-    Uniform Random Number Generator,
C-    giving the same sequence as the IBM and VAX version

      EQUIVALENCE (AMAN,MANT)
      SAVE  MCGN
      DATA  MCGN  /12345/

      MCGN = MCGN * 69069
#if defined(CERNLIB_BUGLRSHFT)
      MANT = ishft (MCGN, -8)
#endif
#if !defined(CERNLIB_BUGLRSHFT)
      MANT = lrshft (MCGN, 8)
#endif
      IF (MANT.EQ.0)         GO TO 14
      AMAN = MANT
C-    AMAN in the range 1 to 2**24-1
C-  - MANT = MANT - X'0C000000'
      MANT = MANT -  201326592
C-    multiply by 2.**(-24)
      RNDM = AMAN
      RETURN

C--   for zero set RNDM = 2.**(-25)
C- 14 MANT = X'33000000'
   14 MANT =  855638016
      RNDM = AMAN
      RETURN

C--       Integer random number
      ENTRY IRNDM (ISEED)
      MCGN = MCGN * 69069
C-  - IRNDM = IAND (MCGN,X'7FFFFFFF')
      IRNDM = AND (MCGN, 2147483647)
      RETURN

C--       Set the seed
      ENTRY RDMIN (ISEED)
      MCGN = ISEED
      RETURN

C--       Get the seed
      ENTRY RDMOUT (ISEED)
      ISEED = MCGN
      END
