*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:13  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GFMATE(IMAT,NAMATE,A,Z,DENS,RADL,ABSL,UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Return parameters for material IMAT                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GTRACK                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER NAMATE*(*)
      DIMENSION UBUF(1)
C.
C.    ------------------------------------------------------------------
C.
      A=-1.
      IF (JMATE.LE.0)GO TO 99
      IF (IMAT.LE.0)GO TO 99
      IF (IMAT.GT.NMATE)GO TO 99
C
      JMA = LQ(JMATE- IMAT)
      IF (JMA.LE.0)GO TO 99
      CALL UHTOC(IQ(JMA+1),4,NAMATE,MIN(20,LEN(NAMATE)))
      A = Q(JMA + 6)
      Z = Q(JMA + 7)
      DENS = Q(JMA + 8)
      RADL = Q(JMA + 9)
      ABSL = Q(JMA + 10)
      NWBUF = IQ(JMA-1) - 11
      IF(NWBUF.GT.0) CALL UCOPY(Q(JMA+12),UBUF,NWBUF)
C
  99  RETURN
      END
