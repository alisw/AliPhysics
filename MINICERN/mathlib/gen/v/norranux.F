*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_UNIX)||defined(CERNLIB_QMALPH)
      FUNCTION NORRAN ( RN )
C
C     FORTRAN version for Alliant of CERN library routine NORRAN
C     (entry V101) for the generation of standard normal pseudo-
C     random numbers.                             / Ch.Walck 880407
C
C     Calling sequences:
C     R = UNI ( DUMMY )              Continuous uniform r.n. 0 TO 1
C     R = VNI ( DUMMY )              Continuous uniform r.n. -1 to 1
C     CALL NORRAN ( R )              Standard normal r.n.
C     CALL NORRIN ( ISEED1, ISEED2 ) Initialize seeds
C     CALL NORRUT ( ISEED1, ISEED2 ) Access seeds
C
      SAVE MCGN, SRGN, TBL
      PARAMETER ( M20=2**20-1, M24=2**24-1, XM24=2.0**(-24)
     +          , XM28=2.0**(-28) )
      DIMENSION TBL(0:326)
      INTEGER R0, R1, R2, SEED1, SEED2, SIGN, SRGN
      EQUIVALENCE (MAN,XMAN)
      DATA TBL/   0.0000,   0.0625, 2*0.1250, 4*0.1875, 5*0.2500
     +        ,   0.5625, 5*0.6250, 3*0.8750,   1.1250,   1.4375
     +        , 5*0.0000, 5*0.0625, 4*0.1250, 2*0.1875,   0.2500
     +        , 5*0.3125, 5*0.3750, 5*0.4375, 5*0.5000, 4*0.5625
     +        , 4*0.6875, 4*0.7500, 4*0.8125,   0.8750, 3*0.9375
     +        , 3*1.0000, 3*1.0625, 2*1.1250, 2*1.1875, 2*1.2500
     +        , 2*1.3125, 2*1.3750,   1.4375,   1.5000,   1.5625
     +        ,   1.6250,   1.6875,   1.7500,   1.8125,10*0.3125
     +        , 7*0.3750, 5*0.4375, 2*0.5000, 9*0.6875, 5*0.7500
     +        ,   0.8125,10*0.9375, 7*1.0000, 3*1.0625,12*1.1875
     +        , 9*1.2500, 5*1.3125, 2*1.3750,13*1.5000,10*1.5625
     +        , 7*1.6250, 5*1.6875, 2*1.7500,15*1.8750,13*1.9375
     +        ,12*2.0000,10*2.0625, 9*2.1250, 8*2.1875, 7*2.2500
     +        , 6*2.3125, 5*2.3750, 4*2.4375, 3*2.5000, 3*2.5625
     +        , 2*2.6250, 2*2.6875 /
      DATA MCGN/12345/, SRGN/1073/
C
      R0   = IEOR ( ISHFT(SRGN,-15), SRGN )
      R1   = ISHFT(R0,17)
      SRGN = IEOR(R0,R1)
      MCGN = 69069 * MCGN
      R2   = IEOR(SRGN,MCGN)
      R0   = ISHFT(R2,-24)
      IF ( R0 .GE. 104 ) GO TO 2
      MAN  = IAND(R2,M24)
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      RN   = XM28 * XMAN + TBL(R0)
      RETURN
C
    2 IF ( R0 .GE. 208 ) GO TO 3
      R0   = R0 - 104
      MAN  = IAND(R2,M24)
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      RN   = - ( XM28 * XMAN + TBL(R0) )
      RETURN
C
    3 R0   = ISHFT(R2,-20)
      IF ( R0 .GE. 3631 ) GO TO 4
      R0   = R0 - 3304
      MAN  = ISHFT ( IAND(R2,M20), 4 )
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      RN   = XM28 * XMAN + TBL(R0)
      RETURN
C
    4 IF ( R0 .GE. 3934 ) GO TO 5
      R0   = R0 - 3607
      MAN  = ISHFT ( IAND(R2,M20), 4 )
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      RN   = - ( XM28 * XMAN + TBL(R0) )
      RETURN
C
    5 RN   = RNORTH ( R2 )
      RETURN
C
C     UNI entry for RNORTH
C
      ENTRY UNI ( DUMMY )
      R0   = IEOR ( ISHFT(SRGN,-15), SRGN )
      R1   = ISHFT ( R0 , 17 )
      SRGN = IEOR ( R0 ,R1 )
      MCGN = 69069 * MCGN
      MAN  = ISHFT ( IEOR ( SRGN, MCGN ) , -8 )
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      UNI = XM24 * XMAN
      RETURN
C
C     Uniform number between -1 and 1 for RNORTH
C
      ENTRY VNI ( DUMMY )
      R0   = IEOR ( ISHFT(SRGN,-15), SRGN )
      R1   = ISHFT ( R0 , 17 )
      SRGN = IEOR ( R0 ,R1 )
      MCGN = 69069 * MCGN
      R2   = IEOR ( SRGN, MCGN )
      MAN  = IAND ( ISHFT ( R2, -7 ), M24 )
      SIGN = ISHFT ( R2, -31 )
      SIGN = 1 - 2 * SIGN
      XMAN = MAN
      IF ( XMAN .EQ. 0.0 ) XMAN = 0.5
      VNI = SIGN * XM24 * XMAN
      RETURN
C
C     Set seeds
C
      ENTRY NORRIN ( SEED1, SEED2 )
      MCGN   = SEED1
      SRGN   = SEED2
      NORRIN = 0
      RETURN
C
C     Get seeds
C
      ENTRY NORRUT ( SEED1, SEED2 )
      SEED1  = MCGN
      SEED2  = SRGN
      NORRUT = 0
      RETURN
      END
#endif
