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
      SUBROUTINE GFNHIT (IUSET, IUDET, NHITS)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GFNHIT (IUSET, IUDET, NHITS*)                         *
C.    *                                                                *
C.    *   Returns the number of hits belonging to the detector IUDET   *
C.    *   of set IUSET                                                 *
C.    *                                                                *
C.    *   IUSET    User set identifier                                 *
C.    *   IUDET    User detector identifier (name of the corresponding *
C.    *            sensivitive volume)                                 *
C.    *   NHITS    returns the total number of hits                    *
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
      NHITS = 0
      IF (JHITS.LE.0)   GO TO 999
      NSET  = IQ(JSET-1)
      CALL GLOOK (IUSET, IQ(JSET+1), NSET, ISET)
      IF (ISET.LE.0)    GO TO 999
*
      JS    = LQ(JSET-ISET)
      JH    = LQ(JHITS-ISET)
      IF (JS.LE.0)      GO TO 999
      IF (JH.LE.0)      GO TO 999
      NDET  = IQ(JS-1)
      CALL GLOOK (IUDET, IQ(JS+1), NDET, IDET)
      IF (IDET.LE.0)    GO TO 999
*
      JD    = LQ(JS-IDET)
      JHD   = LQ(JH-IDET)
      IF (JHD.LE.0)     GO TO 999
      NFIX  = IQ(JD+1) + IQ(JD+3) + 1
      ILAST = IQ(JH+IDET)
*
      IF (ILAST.GT.0) NHITS = ILAST / NFIX
*                                                             END GFNHIT
  999 END
