*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
#include "kerngen/pilot.h"
#if defined(CERNLIB_QMSUN)
#include "sungs/iumode.F"
#elif defined(CERNLIB_QMVAX)
#include "vaxgs/iumode.F"
#else
      FUNCTION IUMODE (WORD)
C
C CERN PROGLIB# M506    IUMODE          .VERSION KERNFOR  4.23  891215
C ORIG. 14/04/79 JZ
C
C     RETURN ZERO IF WORD CONTAINS AN INTEGER
C     PROGRAM WRITTEN 27-APR-79, JZ
C
      INTEGER   WORD
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
C
      J = JBYT (WORD,JEX,NEX)
      IF (J.EQ.NEG)  J=0
      IUMODE= J
      RETURN
      END
#endif
