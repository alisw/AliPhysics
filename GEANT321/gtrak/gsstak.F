*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GSSTAK (IFLAG)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GSSTAK (IFLAG)                                        *
C.    *                                                                *
C.    *   Stores in auxiliary stack JSTAK the particle currently       *
C.    *    described in common /GCKINE/.                               *
C.    *                                                                *
C.    *   On request, creates also an entry in structure JKINE :       *
C.    *    IFLAG =                                                     *
C.    *     0 : No entry in JKINE structure required (user)            *
C.    *     1 : New entry in JVERTX / JKINE structures required (user) *
C.    *    <0 : New entry in JKINE structure at vertex -IFLAG (user)   *
C.    *     2 : Entry in JKINE structure exists already (from GTREVE)  *
C.    *                                                                *
C.    *   Called by : GSKING, GTREVE                                   *
C.    *   Author    : S.Banerjee, F.Bruyant                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmzfo.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
*
      COMMON/VTXKIN/NVTX,ITR
      DIMENSION UBUF(1)
      DATA UBUF/0./
C.
C.    ------------------------------------------------------------------
*
      IF (IPART.LE.0.OR.IPART.GT.NPART) THEN
         PRINT *, ' GSSTAK - Unknown particle code, skip track ', IPART
         GO TO 999
      ENDIF
*
* *** Give control to user for track selection
*
#if !defined(CERNLIB_USRJMP)
      CALL GUSKIP(ISKIP)
#endif
#if defined(CERNLIB_USRJMP)
      CALL JUMPT1(JUSKIP,ISKIP)
#endif
      IF (ISKIP.NE.0) GO TO 999
*
* *** Check if an entry in JKINE structure is required
*
      IF (IFLAG.EQ.1) THEN
         CALL GSVERT (VERT, ITRA, 0, UBUF, 0, NVTX)
         CALL GSKINE (PVERT, IPART, NVTX, UBUF, 0, ITR)
      ELSE IF (IFLAG.LT.0) THEN
         NVTX = -IFLAG
         CALL GSKINE (PVERT, IPART, NVTX, UBUF, 0, ITR)
      ELSE
         IF (IFLAG.EQ.0) THEN
*          Store -ITRA in stack for a track without entry in JKINE
            ITR = -ITRA
         ELSE
            ITR = ITRA
         ENDIF
      ENDIF
*
* *** Store information in stack
*
      IF (JSTAK.EQ.0) THEN
         NDBOOK = NTSTKP*NWSTAK +3
         NDPUSH = NTSTKS*NWSTAK
         CALL MZBOOK (IXCONS,JSTAK,JSTAK,1,'STAK', 0,0,NDBOOK, IOSTAK,3)
         IQ(JSTAK+2) = NTSTKP
      ELSE IF (IQ(JSTAK+1).EQ.IQ(JSTAK+2)) THEN
         CALL MZPUSH (IXCONS, JSTAK, 0, NDPUSH, 'I')
         IQ(JSTAK+2) = IQ(JSTAK+2) +NTSTKS
      ENDIF
*
      JST = JSTAK +IQ(JSTAK+1)*NWSTAK +3
      IQ(JSTAK+1) = IQ(JSTAK+1) +1
      IF (IQ(JSTAK+3).EQ.0) IQ(JSTAK+3) = IQ(JSTAK+1)
      IF (IQ(JSTAK+1).GT.NSTMAX)  NSTMAX = IQ(JSTAK+1)
*
      IQ(JST+1)   = ITR
      IQ(JST+2)   = IPART
      IQ(JST+3)   = 0
      DO 90 I = 1,3
         Q(JST+3+I) = VERT(I)
         Q(JST+6+I) = PVERT(I)
   90 CONTINUE
      Q(JST+10) = TOFG
      Q(JST+11) = SAFETY
      Q(JST+12) = UPWGHT
*
      NALIVE = NALIVE +1
*                                                             END GSSTAK
  999 END
