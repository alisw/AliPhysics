*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1999/05/18 15:55:21  fca
* AliRoot sources
*
* Revision 1.1.1.1  1995/10/24 10:21:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 07/10/94  18.07.13  by  S.Giani
*-- Author :
      SUBROUTINE GTREVE_ROOT
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
#include "geant321/sckine.inc"
      REAL UBUF(2)
      EQUIVALENCE (UBUF(1),WS(1))
      LOGICAL   BTEST
      DIMENSION PMOM(3),VPOS(3)
C.
C.    ------------------------------------------------------------------
      NTMSTO = 0
      NSTMAX = 0
      NALIVE = 0
*         Kick start the creation of the vertex
      VPOS(1)=0
      VPOS(2)=0
      VPOS(3)=0
      PMOM(1)=0
      PMOM(2)=0
      PMOM(3)=0
      IPART=1
      CALL GSVERT(VPOS,0,0,UBUF,0,NVTX)
      CALL GSKINE(PMOM,IPART,NVTX,UBUF,0,NT)
*
      MTRACK=-999
 10   MTROLD=MTRACK
      CALL RXGTRAK(MTRACK,IPART,PMOM,E,VPOS,TTOF)
      IF(MTROLD.LT.0) THEN
         MPRIMA=MTRACK
      ENDIF
      IF(MTRACK.LE.MPRIMA) THEN
         IF(ISWIT(4).GT.0.AND.MTRACK.GT.0) THEN
            IF(ISWIT(4).EQ.1.OR.MOD(MTRACK,ISWIT(4)).EQ.0) THEN
               WRITE(CHMAIL,10200) MTRACK
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
         IF(MTROLD.GT.0) THEN
C --- Output root hits tree only for each primary MTRACK
            CALL RXOUTH
         ENDIF
      ENDIF
      IF(MTRACK.LE.0) GOTO 999
* Set the vertex
      JV=LQ(JVERTX-1)
      Q(JV + 1) = VPOS(1)
      Q(JV + 2) = VPOS(2)
      Q(JV + 3) = VPOS(3)
      Q(JV + 4) = TTOF
      Q(JV + 5) = 0
      Q(JV + 6) = 0
* Set the track
      JK=LQ(JKINE-1)
      Q(JK + 1) = PMOM(1)
      Q(JK + 2) = PMOM(2)
      Q(JK + 3) = PMOM(3)
      Q(JK + 4) = E
      Q(JK + 5) = IPART
      Q(JK + 6) = 1
* Now transport
C      CALL GPVERT(0)
C      CALL GPKINE(0)
* Normal Gtreve_root code
      NV = NVERTX
      DO 40  IV = 1,NV
* ***   For each vertex in turn ..
         JV = LQ(JVERTX-IV)
         NT = Q(JV+7)
         IF (NT.LE.0) GO TO 40
         TOFG   = Q(JV+4)
         SAFETY = 0.
         IF (NJTMAX.GT.0) THEN
            CALL GMEDIA (Q(JV+1), NUMED)
            IF (NUMED.EQ.0) THEN
               WRITE (CHMAIL, 10000) (Q(JV+I), I=1,3)
               CALL GMAIL (0, 0)
               GO TO 40
            ENDIF
            CALL GLSKLT
         ENDIF
*  **   Loop over tracks attached to current vertex
         DO 20  IT = 1,NT
            JV   = LQ(JVERTX-IV)
            ITRA = Q(JV+7+IT)
            IF (BTEST(IQ(LQ(JKINE-ITRA)),0)) GO TO 20
            CALL GFKINE (ITRA, VERT, PVERT, IPART, IVERT, UBUF, NWBUF)
            IF (IVERT.NE.IV) THEN
               WRITE (CHMAIL, 10100) IV, IVERT
               CALL GMAIL (0, 0)
               GO TO 999
            ENDIF
*   *      Store current track parameters in stack JSTAK
            CALL GSSTAK (2)
   20    CONTINUE
*  **   Start tracking phase
   30    IF (NALIVE.NE.0) THEN
            NALIVE = NALIVE -1
*   *      Pick-up next track in stack JSTAK, if any
*   *         Initialize tracking parameters
            CALL GLTRAC
            IF (NUMED.EQ.0) GO TO 30
*   *       Resume tracking
            CALL GUTRAK
            IF (IEOTRI.NE.0) GO TO 999
            GO TO 30
         ENDIF
*
   40 CONTINUE
      GOTO 10
*
10000 FORMAT (' GTREVE_ROOT : Vertex outside setup, XYZ=',3G12.4)
10100 FORMAT (' GTREVE_ROOT : Abnormal track/vertex connection',2I8)
10200 FORMAT (' GTREVE_ROOT : Transporting primary track No ',I10)
*                                                             END GTREVE_ROOT
  999 END
 
