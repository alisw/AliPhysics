*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 17/10/94  10.37.07  by  S.Giani
*-- Author :
      SUBROUTINE GFMOTH(NAME,MOTHER,KONLY)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Fetch MOTHER and ONLY-flag of volume NAME                *
C.    *                                                                *
C.    *    ==>Called by : no more used                                 *
C.    *         Author  P.Zanarini  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
      COMMON/FMOTH/INGLOB,IVOMGL
      CHARACTER*4 NAME,MOTHER
C.
C.    ------------------------------------------------------------------
C.
      MOTHER=' '
      CALL GLOOK(NAME,IQ(JVOLUM+1),NVOLUM,IVOSON)
      IF(IVOSON.LE.0)GO TO 999
      DO 20 IVOMOT=1,NVOLUM
        IVOMGL=IVOMOT
        JVO=LQ(JVOLUM-IVOMOT)
        NIN=Q(JVO+3)
        MNIN=NIN
        IF (NIN.LT.0) MNIN=1
        DO 10 IN=1,MNIN
          INGLOB=IN
          JIN = LQ(JVO-IN)
          IVO=Q(JIN+2)
          IF (IVO.EQ.IVOSON) GO TO 30
   10   CONTINUE
   20 CONTINUE
   30 CALL UHTOC(IQ(JVOLUM+IVOMOT),4,MOTHER,4)
      KONLY=1
      IF (NIN.GT.0) KONLY=Q(JIN+8)
 999  RETURN
      END
