*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFPATH (ISET, IDET, NUMBV, NLEV, LNAM, LNUM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Return the lists of NLEV volume names (LNAM) and numbers  *
C.    *      (LNUM) which identify the path through the JVOLUM data    *
C.    *      structure for the volume corresponding to the detector    *
C.    *      at position IDET in set at position ISET identified by    *
C.    *      the list of node identifiers given in NUMBV.              *
C.    *                                                                *
C.    *      NLEV is expected to be greater than 1                     *
C.    *      and no checks are performed on the validity of ISET/IDET  *
C.    *      In case of user error NLEV returns with the value 0.      *
C.    *                                                                *
C.    *   ==> Called by : <USER>                                       *
C.    *         Author  F.Bruyant  **********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      INTEGER  LNAM(*), LNUM(*), NUMBV(*)
C.
C.    -----------------------------------------------------------------
C.
      JS = LQ(JSET-ISET)
      JD = LQ(JS-IDET)
      NV = IQ(JD+2)
      IPJD = JD +10 +2*NV
      LNAM(1) = IQ(IPJD+1)
      LNUM(1) = 1
      NLEV = IQ(IPJD+2)
C
      NSOL = IQ(JD+9)
      IF (NSOL.EQ.1) THEN
C
C       Usual case
C
        I = 0
        DO 10 N = 2,NLEV
          IPJD = IPJD +2
          LNAM(N) = IQ(IPJD+1)
          LNUM(N) = 1
          IF (IQ(IPJD+2).LE.1) GO TO 10
          I = I +1
          LNUM(N) = NUMBV(I)
   10   CONTINUE
C
      ELSE IF (NSOL.GT.1) THEN
C
C       Case with multiple path
C
        IPSTO = IPJD
        DO 90 IS = 1,NSOL
          IPJDD = JD +8
          DO 30 I = 1,NV
            IPJDD = IPJDD +2
            IF (NUMBV(I).EQ.0)  GO TO 30
            IPJD = IPSTO
            DO 20 N = 2,NLEV
              IPJD = IPJD +2
              IF (IQ(IPJD+1).EQ.IQ(IPJDD+1))  GO TO 30
   20       CONTINUE
            GO TO 81
   30     CONTINUE
C
C         Fill LNAM,LNUM
C
          IPJD = IPSTO
          DO 40 N = 2,NLEV
            IPJD = IPJD +2
            LNAM(N) = IQ(IPJD+1)
            LNUM(N) = 1
   40     CONTINUE
          IPJDD = JD +8
          DO 60 I = 1,NV
            IPJDD = IPJDD +2
            IF (NUMBV(I).EQ.0)  GO TO 60
            IPJD = IPSTO
            DO 50 N = 2,NLEV
              IPJD = IPJD +2
              IF (IQ(IPJD+1).NE.IQ(IPJDD+1))  GO TO 50
              IF (NUMBV(I).GT.IQ(IPJD+2)) GO TO 991
              LNUM(N) = NUMBV(I)
              GO TO 60
   50       CONTINUE
   60     CONTINUE
          GO TO 999
C
   81     IF (IS.EQ.NSOL)  GO TO 991
          IPSTO = IPSTO +2*NLEV
          NLEV = IQ(IPSTO+2)
C
   90   CONTINUE
C
      ELSE
C
C       User error
C
        GO TO 991
C
      ENDIF
      GO TO 999
C
  991 NLEV = 0
      WRITE (CHMAIL, 1000) IQ(JD+9)
      CALL GMAIL(0,0)
C
 1000 FORMAT (' ***** GFPATH USER ERROR, IQ(JD+9)=',I2)
C
  999 RETURN
      END
