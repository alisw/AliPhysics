*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDNSON(NAME,NSON,NDIV)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Return in NSON the number of sons of object NAME         *
C.    *       and in NDIV its number of divisions, always positive     *
C.    *       (NDIV=0 means no division)                               *
C.    *                                                                *
C.    *    ==>Called by : GDPRTR, GDSPEC, GDTR1                        *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 NAME
C.
C.    ------------------------------------------------------------------
C.
      NSON=0
      NDIV=0
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVO)
      IF (IVO.LE.0) GO TO 999
      JVO=LQ(JVOLUM-IVO)
      NIN=Q(JVO+3)
      IF (NIN.EQ.0) GO TO 999
      IF (NIN.GT.0) GO TO 10
C
C             Case of division
C
      NSON=1
      JIN=LQ(JVO-1)
      NDIV=ABS(Q(JIN+3))
      GO TO 999
   10 NSON=NIN
  999 RETURN
      END
