*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:09  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GFNDIG (IUSET, IUDET, NDIGI)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GFNDIG (IUSET, IUDET, NDIGI*)                         *
C.    *                                                                *
C.    *   Returns the number of digits belonging to the detector IUDET *
C.    *   of set IUSET                                                 *
C.    *                                                                *
C.    *   IUSET    User set identifier                                 *
C.    *   IUDET    User detector identifier (name of the corresponding *
C.    *            sensivitive volume)                                 *
C.    *   NDIGI    returns the total number of digits                  *
C.    *                                                                *
C.    *   Called by : <USER>                                           *
C.    *   Author    : S.Banerjee                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      CHARACTER*4     IUSET, IUDET
C.
C.    ------------------------------------------------------------------
*
* *** Find if selected set, detector exists
*
      NDIGI = 0
      IF (JDIGI.LE.0)   GO TO 999
      NSET  = IQ(JSET-1)
      CALL GLOOK (IUSET, IQ(JSET+1), NSET, ISET)
      IF (ISET.LE.0)    GO TO 999
*
      JS    = LQ(JSET-ISET)
      JDI   = LQ(JDIGI-ISET)
      IF (JS.LE.0)      GO TO 999
      IF (JDI.LE.0)     GO TO 999
      NDET  = IQ(JS-1)
      CALL GLOOK (IUDET, IQ(JS+1), NDET, IDET)
      IF (IDET.LE.0)    GO TO 999
*
      JDID  = LQ(JDI-IDET)
      IF (JDID.LE.0)    GO TO 999
      ILAST = IQ(JDI+IDET)
*
      IF (ILAST.NE.0) NDIGI = IQ(JDID+ILAST)
*                                                             END GFNDIG
  999 END
