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
      SUBROUTINE GSKING (IGK)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    SUBR. GSKING (IGK)                                          *
C.    *                                                                *
C.    *   Stores in stack JSTAK either the IGKth track of /GCKING/,    *
C.    *    or the NGKINE tracks when IGK is 0.                         *
C.    *                                                                *
C.    *   Called by : 'User'                                           *
C.    *   Authors   : R.Brun, F.Bruyant                                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
#include "geant321/gcbank.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcnum.inc"
#include "geant321/gctrak.inc"
C.
      COMMON/VTXKIN/NVTX,ITR
      REAL PSTO(3), VSTO(3)
C.
C.    ------------------------------------------------------------------
*
      IF (NGKINE.LE.0) GO TO 999
      IF (IGK.EQ.0) THEN
         N1 = 1
         N2 = NGKINE
      ELSE
         IF (IGK.LT.0.OR.IGK.GT.NGKINE) THEN
            WRITE(CHMAIL,10000) IGK, NGKINE
10000 FORMAT(' GSKING - Abnormal request ',I4,1X,I4)
            CALL GMAIL(0,0)
            GO TO 999
         ENDIF
         N1 = IGK
         N2 = IGK
      ENDIF
*
* *** Save original information from /GCKINE/ and /GCTRAK/
*
      IPASTO = IPART
      TOFSTO = TOFG
      DO 10 I = 1,3
         VSTO(I) = VERT(I)
         PSTO(I) = PVERT(I)
   10 CONTINUE
*
* *** Store required tracks in stack JSTAK
*
      IVSTO = 0
      DO 30 N = N1,N2
         IF (IFLGK(N).LT.0) GO TO 30
         TOFG  = TOFSTO +TOFD(N)
         IPART = GKIN(5,N)
         DO 20 I = 1,3
            VERT(I)  = GPOS(I,N)
            PVERT(I) = GKIN(I,N)
   20    CONTINUE
         IF (IFLGK(N).EQ.1) THEN
            IF (IVSTO.EQ.0) THEN
               IVSTO = 1
               IFLAG = 1
            ELSE
               IFLAG = -NVERTX
            ENDIF
         ELSE
            IFLAG = -IFLGK(N)
         ENDIF
         CALL GSSTAK (IFLAG)
         IF (IFLAG.NE.0) THEN
            IFLGK(N) = ITR
            IFLGK(NGKINE+1) = NVTX
         ENDIF
   30 CONTINUE
*
* *** Restore original information in /GCKINE/ and /GCTRAK/
*
      IPART = IPASTO
      TOFG  = TOFSTO
      DO 40 I = 1,3
         VERT(I)  = VSTO(I)
         PVERT(I) = PSTO(I)
   40 CONTINUE
*                                                             END GSKING
  999 END
