*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDSSUR(P1,P2,P3,P4)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Store space points P1,P2,P3,P4 in vector SURF            *
C.    *       starting from SURF(ISURF)                                *
C.    *                                                                *
C.    *    ==>Called by : GDRWSC                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcvolu.inc"
      DIMENSION SURF(2)
      EQUIVALENCE (WS(1),SURF(1))
      DIMENSION P1(1),P2(1),P3(1),P4(1),RP1(3),RP2(3),RP3(3),RP4(3)
C
C             MPOINT is maximum number of points that can  be stored;
C             at present MPOINT=(2000-12), because SURF is equivalenced
C             from WS(1) to WS(2000)
C
      SAVE MPOINT
      DATA MPOINT/1988/
C.
C.    ------------------------------------------------------------------
C.
C             Rotation matrix
C
      CALL GINROT(P1(1),GRMAT(1,NLEVEL),RP1(1))
      CALL GINROT(P2(1),GRMAT(1,NLEVEL),RP2(1))
      CALL GINROT(P3(1),GRMAT(1,NLEVEL),RP3(1))
      CALL GINROT(P4(1),GRMAT(1,NLEVEL),RP4(1))
C
C             Translation matrix
C
      DO 10 I=1,3
         RP1(I)=RP1(I)+GTRAN(I,NLEVEL)
         RP2(I)=RP2(I)+GTRAN(I,NLEVEL)
         RP3(I)=RP3(I)+GTRAN(I,NLEVEL)
         RP4(I)=RP4(I)+GTRAN(I,NLEVEL)
   10 CONTINUE
C
C             Store space points
C
      SURF(ISURF)=RP1(1)
      SURF(ISURF+1)=RP1(2)
      SURF(ISURF+2)=RP1(3)
      SURF(ISURF+3)=RP2(1)
      SURF(ISURF+4)=RP2(2)
      SURF(ISURF+5)=RP2(3)
      SURF(ISURF+6)=RP3(1)
      SURF(ISURF+7)=RP3(2)
      SURF(ISURF+8)=RP3(3)
      SURF(ISURF+9)=RP4(1)
      SURF(ISURF+10)=RP4(2)
      SURF(ISURF+11)=RP4(3)
      ISURF=ISURF+12
      NSURF=NSURF+1
C
C             If there is not enough space perform a cut
C             to empty the surface buffer SURF
C
      IF (ISURF.GE.MPOINT) CALL GDCUT
C
      END
