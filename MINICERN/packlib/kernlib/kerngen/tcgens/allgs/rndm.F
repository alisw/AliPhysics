*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:50:23  mclareni
* Kernlib
*
*
      FUNCTION RNDM (ISEED)
C
C CERN PROGLIB# V104    RNDM            .VERSION KERNVMI  1.08  930527
C ORIG.  2/02/89  M.K.Storr from IBM version
C
C-    Uniform Random Number Generator,
C-    giving the same sequence as the IBM and VAX version

      REAL         IRNDM
      EQUIVALENCE (AMAN,MANT)
      SAVE  MCGN
      DATA  MCGN  /12345/
      DATA MASK1  /x'0C000000'/, MASK2/x'33000000'/

      MCGN = MCGN * 69069
      MANT = ishft (MCGN,-8)
      IF (MANT.EQ.0)         GO TO 14
      AMAN = MANT
C-    AMAN in the range 1 to 2**24-1
      MANT = MANT - MASK1
C-    multiply by 2.**(-24)
      RNDM = AMAN
      RETURN

C--   for zero set RNDM = 2.**(-25)
   14 MANT = MASK2
      RNDM = AMAN
      RETURN

C--       Integer random number
      ENTRY IRNDM (ISEED)
      MCGN = MCGN * 69069
      MANT = ISHFT (MCGN,-1)
      IRNDM = AMAN
      RETURN

C--       Set the seed
      ENTRY RDMIN (ISEED)
      MCGN = ISEED
      RNDM  = 0.
      RETURN

C--       Get the seed
      ENTRY RDMOUT (ISEED)
      ISEED = MCGN
      RNDM  = 0.
      END
