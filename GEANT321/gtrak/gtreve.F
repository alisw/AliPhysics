*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 07/10/94  18.07.13  by  S.Giani
*-- Author :
      SUBROUTINE GTREVE
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GTREVE                                                *
C.    *                                                                *
C.    *   Controls tracking of all particles belonging to the current  *
C.    *    event.                                                      *
C.    *                                                                *
C.    *   Called by : GUTREV, called by GTRIG                          *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#include "geant321/gcvdma.inc"
      REAL UBUF(2)
      EQUIVALENCE (UBUF(1),WS(1))
      LOGICAL   BTEST
C.
C.    ------------------------------------------------------------------
      NTMSTO = 0
      NSTMAX = 0
      NALIVE = 0
*
      NV = NVERTX
      DO 290 IV = 1,NV
*
* ***   For each vertex in turn ..
*
         JV = LQ(JVERTX-IV)
         NT = Q(JV+7)
         IF (NT.LE.0) GO TO 290
         TOFG   = Q(JV+4)
         SAFETY = 0.
         IF (NJTMAX.GT.0) THEN
            CALL GMEDIA (Q(JV+1), NUMED)
            IF (NUMED.EQ.0) THEN
               WRITE (CHMAIL, 1001) (Q(JV+I), I=1,3)
               CALL GMAIL (0, 0)
               GO TO 290
            ENDIF
            CALL GLSKLT
         ENDIF
*
*  **   Loop over tracks attached to current vertex
*
         DO 190 IT = 1,NT
            JV   = LQ(JVERTX-IV)
            ITRA = Q(JV+7+IT)
            IF (BTEST(IQ(LQ(JKINE-ITRA)),0)) GO TO 190
            CALL GFKINE (ITRA, VERT, PVERT, IPART, IVERT, UBUF, NWBUF)
            IF (IVERT.NE.IV) THEN
               WRITE (CHMAIL, 1002) IV, IVERT
               CALL GMAIL (0, 0)
               GO TO 999
            ENDIF
*
*   *      Store current track parameters in stack JSTAK
*
            CALL GSSTAK (2)
  190    CONTINUE
*
*  **   Start tracking phase
*
  210    IF (NALIVE.NE.0) THEN
            NALIVE = NALIVE -1
*
*   *      Pick-up next track in stack JSTAK, if any
*
            IF (IQ(JSTAK+1).GT.0) THEN
*
*   *         Initialize tracking parameters
*
               CALL GLTRAC
               IF (NUMED.EQ.0) GO TO 210
            ELSE
*
*   *         otherwise, select next track segment from stack JTRACK
*
               CALL GFTRAC
*
            ENDIF
*
*   *       Resume tracking
*
            IF(RAYTRA.EQ.1.)THEN
#if defined(CERNLIB_HIGZ)
              CALL GDRACK
#endif
            ELSE
#if !defined(CERNLIB_USRJMP)
              CALL GUTRAK
#endif
#if defined(CERNLIB_USRJMP)
              CALL JUMPT0(JUTRAK)
#endif
            ENDIF
            IF (IEOTRI.NE.0) GO TO 999
            GO TO 210
         ENDIF
*
  290 CONTINUE
*
 1001 FORMAT (' GTREVE : Vertex outside setup, XYZ=',3G12.4)
 1002 FORMAT (' GTREVE : Abnormal track/vertex connection',2I8)
*                                                             END GTREVE
  999 END
 
