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
C CERN PROGLIB# V104    RNDM            .VERSION KERNVAX  2.39  930112
C ORIG. 22-MAR-88 from IBM version of G.Marsaglia et al., Montreal
C
C     !!!   Version valid only for VAX floating point format  !!!
C     Careful: compile with /CHECK=NOOVERFLOW
C
C-    Uniform Random Number Generator,
C-    giving the same sequence as the IBM and VAX m-language version

      EQUIVALENCE (AMAN,MANT)
      PARAMETER    (MSK1 = 3072)
C     PARAMETER    (MSK1 = '00000C00'X)
      PARAMETER    (MSK2 = 13312)
C     PARAMETER    (MSK2 = '00003400'X)
      SAVE  MCGN
      DATA  MCGN  /12345/


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
      RETURN

C--       Set the seed
      ENTRY RDMIN (ISEED)
      MCGN = ISEED
      RETURN

C--       Get the seed
      ENTRY RDMOUT (ISEED)
      ISEED = MCGN
      END
