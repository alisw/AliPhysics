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
      SUBROUTINE GDSON(N,NAME,ISON)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Returns in ISON the N-th son of object NAME              *
C.    *                                                                *
C.    *    ==>Called by : GDPRTR, GDSPEC, GDTR1                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 NAME,ISON
C.
C.    ------------------------------------------------------------------
C.
      ISON=' '
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF (IVO.LE.0) GO TO 999
      JVO=LQ(JVOLUM-IVO)
      NIN=Q(JVO+3)
      IF (NIN.EQ.0) GO TO 999
      IF (NIN.GT.0) GO TO 10
      IF (N.NE.1) GO TO 999
      JIN=LQ(JVO-1)
      GO TO 20
   10 JIN=LQ(JVO-N)
   20 IVOL=Q(JIN+2)
      CALL UHTOC(IQ(JVOLUM+IVOL),4,ISON,4)
  999 RETURN
      END
