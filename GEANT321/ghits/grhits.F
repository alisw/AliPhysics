*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:11  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GRHITS (IUSET, IUDET, NTRA, ITRA)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GRHITS (IUSET, IUDET, NTRA, ITRA)                     *
C.    *                                                                *
C.    *   Removes the hits produced by the tracks ITRA(1...NTRA) in    *
C.    *   the generic volume name IUDET belonging to the set IUSET.    *
C.    *                                                                *
C.    *   IUSET    User set identifier                                 *
C.    *   IUDET    User detector identifier (name of the corresponding *
C.    *            sensivitive volume)                                 *
C.    *   NTRA     Number of tracks whose hits are to be removed       *
C.    *   ITRA     Track indices whose hits are to be removed          *
C.    *                                                                *
C.    *   Called by : <USER>                                           *
C.    *   Author    : S.Banerjee                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      COMMON /GC1HIT/ LOC(2), JD, JDH, JH, JHD, JHDN, JS
      DIMENSION       ITRA(*)
      CHARACTER*(*)   IUSET, IUDET
      SAVE JHDNN
      DATA JHDNN/0/
*
*     ------------------------------------------------------------------
*
      IF (NTRA.LE.0)                  GO TO 999
      IF (JSET.LE.0)                  GO TO 999
      IF (JHITS.LE.0)                 GO TO 999
*
* *** Find the selected set(s)
*
      NSET  = IQ(JSET-1)
      IF (IUSET(1:1).EQ.'*') THEN
         NS1   = 1
         NS2   = NSET
      ELSE
         CALL GLOOK (IUSET, IQ(JSET+1), NSET, ISET)
         IF (ISET.LE.0)               GO TO 999
         NS1   = ISET
         NS2   = ISET
      ENDIF
      CALL MZLINT (IXSTOR, '/GC1HIT/', LOC, JD, JS)
*
* *** Loop over selected sets
*
      DO 30  ISET = NS1, NS2
         JS    = LQ(JSET-ISET)
         JH    = LQ(JHITS-ISET)
         IF (JS.LE.0.OR.JH.LE.0)      GO TO 30
         NDET  = IQ(JS-1)
*
*  **    Find the selected detector(s)
*
         IF (IUDET(1:1).EQ.'*') THEN
            ND1   = 1
            ND2   = NDET
         ELSE
            CALL GLOOK (IUDET, IQ(JS+1), NDET, IDET)
            IF (IDET.EQ.0)            GO TO 30
            ND1   = IDET
            ND2   = IDET
         ENDIF
*
*  **    Loop over selected detectors
*
         DO 20 IDET = ND1, ND2
            JD    = LQ(JS-IDET)
            JHD   = LQ(JH-IDET)
            IF (JD.LE.0)              GO TO 20
            IF (JHD.LE.0)             GO TO 20
            JDH   = LQ(JD-1)
            IF (JDH.LE.0)             GO TO 20
            ILAST = IQ(JH+IDET)
            IF (ILAST.EQ.0)           GO TO 20
            NW    = IQ(JD+1) + IQ(JD+3) + 1
*
*  **       Shunt the original bank and lift a new SJDH bank
*
            CALL ZSHUNT (IXDIV, JHD, JHDNN, 2, 0)
            CALL MZBOOK (IXDIV, JHDN, JH, -IDET, 'SJHD', 0, 0, ILAST,
     +                   1, -1)
            IQ(JHDN-5) = IQ(JHD-5)
*
*  **       Copy the relevant part
*
            LAST  = 0
            DO 10 I = 1, ILAST, NW
               II    = IUCOMP (IQ(JHD+I), ITRA, NTRA)
               IF (II.LE.0) THEN
                  CALL UCOPY (IQ(JHD+I), IQ(JHDN+LAST+1), NW)
                  LAST  = LAST + NW
               ENDIF
   10       CONTINUE
*
*  **       Drop the old bank
*
            IF (LAST.LT.ILAST) THEN
               CALL VZERO (IQ(JHDN+LAST+1), ILAST-LAST)
            ENDIF
            IQ(JH+IDET) = LAST
            CALL MZDROP (IXDIV, JHD, ' ')
   20    CONTINUE
   30 CONTINUE
*
  100 LOC(1) = 0
*                                                             END GRHITS
  999 END
