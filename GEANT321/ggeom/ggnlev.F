*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GGNLEV
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GGNLEV                                                *
C.    *                                                                *
C.    *   Scans the JVOLUM structure to find the maximum Number of     *
C.    *   LEVels and maximum number of volumes in each level.top       *
C.    *   Stores corresponding information in bank next to JVOLUM      *
C.    *   mother bank.                                                 *
C.    *                                                                *
C.    *   Called by : GGCLOS                                           *
C.    *   Author    : S.Banerjee                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      PARAMETER (NLVMAX=15)
      INTEGER ILINK(NLVMAX), LJVOM(NLVMAX), NLINK(NLVMAX)
      INTEGER NDVLP(NLVMAX), MAXVO(NLVMAX)
C.
C.    ------------------------------------------------------------------
*
      CALL VZERO (MAXVO, NLVMAX)
      NLEV = 1
      NLVT = NLEV
      LJVOM(NLEV) = LQ(JVOLUM-1)
      ILINK(NLEV) = 1
      MAXVO(NLEV) = 1
*
* *** Loop over volumes going down the tree
*
   10 IF (NLEV.GT.NLVT) NLVT = NLEV
*
*  ** Check if current volume is a source of local development
*
      IF (LQ(LJVOM(NLEV)).GT.0) THEN
         NDVLP(NLEV) = NLEV
      ELSE IF (NLEV.GT.1) THEN
         NDVLP(NLEV) = NDVLP(NLEV-1)
      ELSE
         NDVLP(NLEV) = NLVMAX
      ENDIF
      NIN = Q(LJVOM(NLEV)+3)
      IF (NIN.GT.0) THEN
*
*   *    Current volume has contents defined by position
*
         NLINK(NLEV) = NIN
         IF (NIN.GT.MAXVO(NLEV+1)) MAXVO(NLEV+1) = NIN
         IF (NLINK(NLEV).EQ.0) GO TO 30
         IN   = NLINK(NLEV)
         NLINK(NLEV) = IN - 1
         JIN  = LQ(LJVOM(NLEV)-IN)
         IVOT = Q(JIN+2)
         NLEV = NLEV + 1
         LJVOM(NLEV) = LQ(JVOLUM-IVOT)
         ILINK(NLEV) = IN
         GO TO 10
*
      ELSE IF (NIN.LT.0) THEN
*
*   *    Current volume is divided
*
         JDIV = LQ(LJVOM(NLEV)-1)
         IVOT = Q(JDIV+2)
         IF (NLEV.GT.NDVLP(NLEV)) THEN
            JPAR = LQ(LJVOM(NDVLP(NLEV)))
            DO 15 ILEV = NDVLP(NLEV), NLEV-1
               IF (IQ(JPAR+1).EQ.0) THEN
                  JPAR = LQ(JPAR-ILINK(ILEV+1))
                  IF (JPAR.EQ.0) GO TO 20
               ELSE
                  JPAR = LQ(JPAR-1)
               ENDIF
               IF (ILEV.EQ.NLEV-1) THEN
                  NDIV = IQ(JPAR+1)
                  GO TO 25
               ENDIF
   15       CONTINUE
         ENDIF
   20    NDIV = Q(JDIV+3)
         IF (NDIV.LT.0) NDIV = -NDIV
   25    NLINK(NLEV) = 0
         IF (NDIV.GT.MAXVO(NLEV+1)) MAXVO(NLEV+1) = NDIV
         NLEV = NLEV + 1
         LJVOM(NLEV) = LQ(JVOLUM-IVOT)
         ILINK(NLEV) = 1
         GO TO 10
      ENDIF
*
*  ** No more content in this volume
*   * Go up the tree
*
   30 NLEV = NLEV - 1
      IF (NLEV.GT.0) THEN
         IF (NLINK(NLEV).EQ.0) GO TO 30
         IN   = NLINK(NLEV)
         NLINK(NLEV) = IN - 1
         JIN  = LQ(LJVOM(NLEV)-IN)
         IVOT = Q(JIN+2)
         NLEV = NLEV + 1
         LJVOM(NLEV) = LQ(JVOLUM-IVOT)
         ILINK(NLEV) = IN
         GO TO 10
      ENDIF
*
* *** Scan is complete; now store the information in next-of-same-type
*  ** link of JVOLUM
*
      JVOLX = LQ(JVOLUM)
      IF (JVOLX.GT.0) CALL MZDROP (IXCONS, JVOLX, ' ')
      CALL MZBOOK (IXCONS, JVOLX, JVOLUM, 0, 'NLEV', 0, 0, NLVT, 2, 0)
      CALL UCOPY (MAXVO, IQ(JVOLX+1), NLVT)
#if defined(CERNLIB_DEBUGG)
      IF (IDEBUG.NE.0) THEN
         WRITE (CHMAIL, 2001) NLVT
         CALL GMAIL (1, 0)
         IF (NLVT.GT.0) THEN
            WRITE (CHMAIL, 2002) (MAXVO(I), I = 1, NLVT)
            CALL GMAIL (0, 0)
         ENDIF
      ENDIF
 2001 FORMAT (' GGNLEV : Maximum number of levels = ',I4,
     +        ' -- maximum number of volumes in each level')
 2002 FORMAT (5X,20I6)
#endif
*                                                             END GGNLEV
      END
