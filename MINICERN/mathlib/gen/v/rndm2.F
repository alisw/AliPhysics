*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_VAX)
      OPTIONS /CHECK=NOOVERFLOW
      FUNCTION RNDM2 ()
C-----------------------------------------------------------------------
C
C...       FUNCTION RNDM2  ()                               Ch.Walck    870401
C.                                                        slight mod.   940122
C.
C.    Fortran version of RNDM2/IRNDM2/RD2IN/RD2OUT (CERN LIBRARY entry
C.    V107 on GENLIB) for VAX/VMS or AXP/OpenVMS written after the
C.    assembler routine.
C.
C.    Calling sequences:
C.    R = RNDM2(DUMMY)                    Generate continuous uniform r.n.
C.    IR = IRNDM2(DUMMY)                  Generate discrete uniform r.n.
C.    CALL RD2IN ( ISEED1, ISEED2 )       Initialize seeds
C.    CALL RD2OUT ( ISEED1, ISEED2 )      Access seeds
C.
C.    This generator has a sequence of 4.609.432.020.664.188.928 =
C.    4.6*10**18 random numbers. The minimum value returned is 2**-25 =
C.    0.000000030 and the maximum value is 1 - 2**24 = 0.999999940.
C.
C.    The arithmetics involved causes integer overflows which is handled
C.    on VAX/VMS by calling ERRSET once from the main program as follows:
C.       CALL ERRSET(70,.TRUE.,.FALSE.,.FALSE.,.FALSE.,100)
C.    or by compiling this routine with /CHECK=NOOVERFLOW
C.    or, as here, by the VAX FORTRAN OPTIONS statement before and after
C.       the routine.
C.
C.----------------------------------------------------------------------
      PARAMETER (MEXPO='FFFF807F'X,NEXPO=7)
      SAVE MCGN, SRGN
      INTEGER R0, R1, SRGN
      EQUIVALENCE ( MAN, XMAN ), ( IR, R )
      DATA MCGN/12345/, SRGN/1073/
      R0    = IEOR ( ISHFT(SRGN,-15), SRGN )
      R1    = ISHFT ( R0 , 17 )
      SRGN  = IEOR ( R0 ,R1 )
      MCGN  = 69069 * MCGN
      MAN   = ISHFT ( IEOR ( SRGN, MCGN ) , -8 )
      XMAN  = MAN
      NSHFT = IAND ( ISHFT ( MAN, -NEXPO ), 31 )
      NSHFT = NSHFT + 104
      IR    = IOR ( IAND ( MAN, MEXPO ), ISHFT ( NSHFT, NEXPO ) )
      RNDM2 = R
      RETURN
C
      ENTRY IRNDM2 ()
      R0     = IEOR ( ISHFT(SRGN,-15), SRGN )
      R1     = ISHFT ( R0 , 17 )
      SRGN   = IEOR ( R0, R1 )
      MCGN   = 69069 * MCGN
      IRNDM2 = ISHFT ( IEOR ( MCGN, SRGN ) , -1 )
      RETURN
C
      ENTRY RD2IN ( ISEED1, ISEED2 )
      MCGN = ISEED1
      SRGN = ISEED2
      RETURN
C
      ENTRY RD2OUT ( ISEED1, ISEED2 )
      ISEED1 = MCGN
      ISEED2 = SRGN
      RETURN
      END
#endif
#if (defined(CERNLIB_VAX))&&(defined(CERNLIB_FORTRAN))&&(!defined(CERNLIB_F4))&&(defined(CERNLIB_OLD))
        FUNCTION RNDM2 (DUMMY)
C
C       VERY USEFULL FORTRAN VERSION OF RNDM2
C
        INTEGER*4 M(6)
C
        INTEGER*4 MCGN,SRGN,REGB,REGC,REGD
        DATA MCGN,SRGN /12345,1073/
C
        DATA M /'00F00000'X,'000F0000'X,'0000F000'X,
     1          '00000F00'X,'000000F0'X,'0000000F'X/
C
        DATA ICALL /0/
        ICALL= ICALL+1
*       WRITE(6,100)ICALL
C
        REGB= SRGN
        REGC= REGB
        REGC= ISHFT (REGC,-15)
        REGB= IEOR (REGB,REGC)
        REGC= REGB
        REGC= ISHFT (REGC,17)
        REGB= IEOR (REGB,REGC)
        SRGN= REGB
        REGD= MCGN
        REGD= REGD*69069
        MCGN= REGD
*       WRITE(6,110)MCGN,SRGN
        REGD= IEOR (REGB,REGD)
        REGD= ISHFT (REGD,-8)
        REGD= IAND(REGD,'40000000'X)
        R=0.
        DO II=1,6
          IP=(II-6)*4
          J=IAND(M(II),REGD)
          JJ=ISHFT(J,IP)
          R= R+JJ*(16.**-II)
        ENDDO
        RNDM2= R
        RETURN
 100    FORMAT(5X,I3,'-TH CALL TO RNDM2.FOR ')
 110    FORMAT(5X,'MCGN,SRGN = ',2Z9.8)
        END
#endif
#if defined(CERNLIB_UNIX)
      FUNCTION RNDM2 (DUMMY)
C
C CERN PROGLIB# V107    RNDM2           .VERSION KERNALT  1.02  880323
C ORIG. 7-APR-88, Ch. Walck, Stockholm
C
C Mods       Date   Comments
C MARQUINA 90/05/06 Generalize code for UNIX machines
C
C-    Uniform Random Number Generator

C-    Calling sequences:
C-    R  = RNDM2 ()               Continuous uniform r.n. 0 to 1
C-    IR = IRNDM2 ()              Discrete uniform r.n. 0 to 2**31-1
C-    CALL RD2IN  (ISEED1,ISEED2)    Set seeds
C-    CALL RD2OUT (ISEED1,ISEED2)    Get seeds

      EQUIVALENCE (ZERO,IZERO)
      SAVE  MCGN, MCGX
      DATA  MCGN/12345/, MCGX/1073/
#include "v107z0.inc"

C----      Floating random number

#include "v107rn.inc"
      IF (MANT.EQ.0)         GO TO 14
      AMAN = MANT
      RNDM2 = AMAN * 2.**(-24)
      RETURN

C--   for zero set RNDM2 = 2.**(-25)
   14 IZERO = IZ0
      RNDM2 = ZERO
      RETURN

C---       Integer random number

      ENTRY IRNDM2 (DUMMY)
#include "v107ri.inc"
      RETURN

C----      Set seeds

      ENTRY RD2IN (ISEED1, ISEED2)
      MCGN = ISEED1
      MCGX = ISEED2
      RETURN

C----      Get seeds

      ENTRY RD2OUT (ISEED1, ISEED2)
      ISEED1 = MCGN
      ISEED2 = MCGX
      RETURN
      END
#endif
