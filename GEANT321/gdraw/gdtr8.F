*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR8(NAME,MOTH,IONL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Like routine GFMOTH                                      *
C.    *                                                                *
C.    *       NAME = volume name (input)                               *
C.    *       MOTH = mother name of NAME (output)                      *
C.    *       IONL = 'ONLY' flag (1 for 'ONLY', 0 for 'MANY') (output) *
C.    *                                                                *
C.    *    ==>Called by : GDTREE                                       *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      CHARACTER*4 NAME,MOTH
C.
C.    ------------------------------------------------------------------
C.
      IONL=1
      CALL GLOOK(NAME,IQ(JNAM1+1),NUMND1,J)
      JM=IQ(JMOT1+J)
      IF (JM.EQ.0) GO TO 999
      CALL UHTOC(IQ(JNAM1+JM),4,MOTH,4)
      CALL GLOOK(MOTH,IQ(JVOLUM+1),NVOLUM,IVOMOT)
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVOSON)
      JVO=LQ(JVOLUM-IVOMOT)
      NIN=Q(JVO+3)
      MNIN=NIN
      IF (NIN.LT.0) MNIN=1
      DO 10 IN=1,MNIN
         JIN=LQ(JVO-IN)
         IVO=Q(JIN+2)
         IF (IVO.EQ.IVOSON) GO TO 20
   10 CONTINUE
      GO TO 999
   20 CONTINUE
      IF (NIN.GT.0) IONL=Q(JIN+8)
  999 RETURN
      END
 
