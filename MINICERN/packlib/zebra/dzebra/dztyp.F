*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:07  mclareni
* Zebra
*
*
*------------------------------------------------------------
#include "zebra/pilot.h"
#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf1.inc"
#endif
      SUBROUTINE DZTYP
      SAVE ALOW,AUP,BLOW,BUP,ILIM
#include "zebra/mqsys.inc"
#include "zebra/qequ.inc"
#include "zebra/zbcd.inc"
#include "zebra/zbcdk.inc"
#include "zebra/zmach.inc"
#include "zebra/dzc1.inc"

      PARAMETER (IUMOIQ = 0)

      INTEGER IA,IWORD
      REAL     A, WORD
      EQUIVALENCE (WORD,IWORD), (A,IA)

      CHARACTER CHROUT*(*),CHSTAK*6
      PARAMETER (CHROUT = 'DZTYP' )

#if (defined(CERNLIB_DEBUGON))&&(defined(CERNLIB_VFORT))
#include "zebra/debugvf2.inc"
#endif

      DATA ILIM/10000000/,ALOW,AUP/.001,1.0E6/,BLOW,BUP/1.0E-20,1.0E20/



      CHSTAK          = CQSTAK(MCQSIQ:)
      CQSTAK(MCQSIQ:) = CHROUT

      J      = LBASE + JD
      IWORD  = LQ(KQS+J)
      JTYP   = 1
      IF (IWORD.EQ.IQNIL.OR.IWORD.EQ.IQNIL+J)              GO TO 999

      IDLINK = JD+1 - JDFD

      IIUMOD = IUMODE(WORD)

      IF (IDLINK.LE.0)           THEN
          IF (IIUMOD.NE.IUMOIQ)                            GO TO 999
          IF (IWORD.EQ.LNULL)              THEN
              JTYP   = 2
          ELSE
              CALL MZCHLS(NCHEKQ,IWORD)
              JTYP = -1
          ENDIF
                                                           GO TO 999
      ENDIF



      IF (IFLOPT(MPOSZQ).NE.0)                             GO TO 100


      ITYPE  = 0

      IF(IIUMOD.EQ.IUMOIQ)                                 GO TO 300
      IF (IWORD.EQ.IQBLAN)                                 GO TO 200
      IF (WORD.EQ.0.) THEN
                                                           GO TO 400
      ELSEIF (WORD.LT.0.) THEN
          IF (WORD.LT.-BLOW.AND.WORD.GT.-BUP)              GO TO 400
      ELSE
          IF (WORD.GT.BLOW.AND.WORD.LT.BUP)                GO TO 400
      ENDIF
      CALL UBLOW(WORD,IQUEST,NQCHAW)
      DO 90 J=1,NQCHAW
          IF (IZBCD(IQUEST(J)).GE.48)                      GO TO 100
   90 CONTINUE
                                                           GO TO 200


 100  JTYP   = 1
                                                           GO TO 999


  200 JTYP   = 7
                                                           GO TO 999


  300 IF ((IWORD.GE.0.AND.IWORD.LT.ILIM).OR.
     X    (IWORD.LT.0.AND.IWORD.GT.-ILIM)   ) THEN
          JTYP   = 2
      ENDIF
                                                           GO TO 999


  400 A      = ABS(WORD)
      JTYP   = 6
      IF     (A.EQ.0.)                            THEN
          JTYP   = 3
      ELSEIF (A.LT.AUP.AND.A.GT.ALOW)             THEN
          JTYP   = 4
          IF (A.EQ.AINT(A))                           THEN
              JTYP   = 3
          ELSEIF (A.LT.100.)                          THEN
              JTYP   = 5
          ENDIF
      ENDIF

  999 CQSTAK(MCQSIQ:) = CHSTAK
      END
