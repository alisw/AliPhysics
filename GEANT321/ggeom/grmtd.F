*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:54  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GRMTD(XMAS,IROT,XNEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE ROTATES A POINT IN XMAS FROM THE MASTER SYSTEM      *
C.    *    TO A NEW SYSTEM USING MATRIX # IROT.                        *
C.    *    RESULT IS RETURNED IN VECTOR XNEW.                          *
C.    *                                                                *
C.    *    ==>Called by : GNEXT, GTNEXT                                *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION XMAS(3),XNEW(3)
C.
C.    ------------------------------------------------------------------
C.
      IF (IROT.GT.0) THEN
         JR = LQ(JROTM- IROT)
         XNEW(1) = XMAS(1)*Q(JR+1) + XMAS(2)*Q(JR+2) + XMAS(3)*Q(JR+3)
         XNEW(2) = XMAS(1)*Q(JR+4) + XMAS(2)*Q(JR+5) + XMAS(3)*Q(JR+6)
         XNEW(3) = XMAS(1)*Q(JR+7) + XMAS(2)*Q(JR+8) + XMAS(3)*Q(JR+9)
*
      ELSE
         XNEW(1) = XMAS(1)
         XNEW(2) = XMAS(2)
         XNEW(3) = XMAS(3)
*
      ENDIF
      END
