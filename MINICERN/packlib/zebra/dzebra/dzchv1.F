*
* $Id$
*
* $Log$
* Revision 1.2  1996/04/24 17:26:07  mclareni
* Extend the include file cleanup to dzebra, rz and tq, and also add
* dependencies in some cases.
*
* Revision 1.1.1.1  1996/03/06 10:47:06  mclareni
* Zebra
*
*
*-----------------------------------------------------------
#include "zebra/pilot.h"
      SUBROUTINE DZCHV1 (LBEGIN,LEND,IXQUES,ISUM)
#if defined(CERNLIB_QMCRY)
CDIR$ INTEGER=64
#endif
#include "zebra/mqsys.inc"
#include "zebra/zbcdch.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zunit.inc"
#include "zebra/dzc1.inc"
#include "zebra/questparq.inc"
#include "zebra/storparq.inc"
      INTEGER ISUM(*)
      PARAMETER ( NFIELD =  4 )
*       32 BIT MACHINES
#if defined(CERNLIB_B32)
      PARAMETER ( NBITS  =  8 )
*       36 BIT MACHINE
#endif
#if defined(CERNLIB_B36)
      PARAMETER ( NBITS  =  9 )
*       60 BIT MACHINE
#endif
#if defined(CERNLIB_B60)
      PARAMETER ( NBITS  = 15 )
*       64 BIT MACHINE
#endif
#if defined(CERNLIB_B64)
      PARAMETER ( NBITS  = 16 )
#endif
      PARAMETER ( NWMAX  = 2**(NBITS*(NFIELD-1)-1) )
      INTEGER IFIELD(NFIELD)

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZCHV1')

#include "zebra/q_jbyt.inc"

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif




      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      NWTOT  = LEND - LBEGIN + 1
      IF (NWTOT.GT.NWMAX)              THEN
          WRITE(CQINFO,'(2I10)') NWTOT,NWMAX
          CALL DZTEXT(MCHV1Q,CDUMMQ,0)
          IQUEST(1) = 1
                                                           GO TO 999
      ENDIF

      IF (IXQUES.NE.0) THEN
          DO 10 JF=1,NFIELD
   10     IFIELD(JF) = IQUEST(IXQUES+JF)
      ELSE
          DO 20 JF=1,NFIELD
   20     IFIELD(JF) = 0
      ENDIF


      DO 100 JW=LBEGIN,LEND
          DO 100 JF=1,NFIELD
              JFIELD=JBYT(LQ(JW),(JF-1)*NBITS+1,NBITS)
  100     IFIELD (JF) = IFIELD(JF) + JFIELD

      IF (IXQUES.EQ.0) THEN
          JCARRY = 0

          DO 200 JF=1,NFIELD
              IFIELD (JF) = IFIELD(JF) + JCARRY
              JCARRY = IFIELD(JF)/2**NBITS
          IFIELD (JF) = IFIELD(JF) - JCARRY*(2**NBITS)
  200     CALL SBYT(IFIELD(JF),ISUM(1),(JF-1)*NBITS+1,NBITS)
          ISUM(2) = JCARRY
      ELSE
          DO 300 JF=1,4
  300     IQUEST(IXQUES+JF) = IFIELD(JF)
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      RETURN
      END
