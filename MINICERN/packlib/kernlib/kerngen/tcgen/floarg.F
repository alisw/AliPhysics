*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:44  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/floarg.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/floarg.F"
#else
      FUNCTION FLOARG (ARGM)
C
C CERN PROGLIB# M250    FLOARG          .VERSION KERNFOR  4.23  891215
C ORIG. 01/8/82  ADAPTED FROM UFLINT
C
C
      DIMENSION    ARGM(9)
      EQUIVALENCE  (A,IA)
C
#if defined(CERNLIB_QMCV64)
      PARAMETER    (JEX=53, NEX=12, NEG=4095 )
#endif
#if defined(CERNLIB_QMCRY)
      PARAMETER    (JEX=49, NEX=16, NEG=65535)
#endif
#if defined(CERNLIB_QMCDC)
      PARAMETER    (JEX=49, NEX=12, NEG=4095 )
#endif
#if defined(CERNLIB_QMUNI)||defined(CERNLIB_QMUNO)
      PARAMETER    (JEX=28, NEX= 9, NEG= 511 )
#endif
#if defined(CERNLIB_QMPDP)
      PARAMETER    (JEX=28, NEX= 9, NEG= 511 )
#endif
#if defined(CERNLIB_QMDGE)||defined(CERNLIB_QMIBM)
      PARAMETER    (JEX=25, NEX= 8, NEG= 255 )
#endif
#if (defined(CERNLIB_QIEEE))&&(defined(CERNLIB_B32))
      PARAMETER    (JEX=24, NEX= 9, NEG= 511 )
#endif
#if defined(CERNLIB_QMND3)
      PARAMETER    (JEX=23, NEX=10, NEG=1023 )
#endif
C
#include "kerngen/q_jbit.inc"
* Ignoring t=pass
C
C
      A   = ARGM(1)
      JCH = JBYT (IA,JEX,NEX)
      IF (JCH.EQ.0)          GO TO 36
      IF (JCH.EQ.NEG)        GO TO 36
      FLOARG = A
      RETURN
C
   36 FLOARG = IA
      RETURN
      END
#endif
