*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:14  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GFTMED(NUMED,NATMED,NMAT,ISVOL,IFIELD,FIELDM,
     +           TMAXFD,STEMAX,DEEMAX,EPSIL,STMIN,UBUF,NWBUF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Return parameters for tracking medium NUMED              *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GTRACK                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER*20 NATMED
      DIMENSION UBUF(1)
C.
C.    ------------------------------------------------------------------
C.
      NMAT=0
      IF (JTMED.LE.0)GO TO 99
      IF (NUMED.LE.0)GO TO 99
      IF (NUMED.GT.NTMED)GO TO 99
      JTM = LQ(JTMED- NUMED)
      IF (JTM.LE.0)GO TO 99
C
      CALL UHTOC(IQ(JTM+1),4,NATMED,20)
      NMAT = Q(JTM + 6)
      ISVOL = Q(JTM + 7)
      IFIELD = Q(JTM + 8)
      FIELDM = Q(JTM + 9)
      TMAXFD = Q(JTM + 10)
      STEMAX = Q(JTM + 11)
      DEEMAX = Q(JTM + 12)
      EPSIL = Q(JTM + 13)
      STMIN = Q(JTM + 14)
      NWBUF = IQ(JTM-1) - 14
      IF(NWBUF.LE.0) GO TO 99
      DO 20 I=1,NWBUF
  20  UBUF(I)=Q(JTM+14+I)
C
  99  RETURN
      END
