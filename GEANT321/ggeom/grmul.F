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
      SUBROUTINE GRMUL(RMAT,IROT,RMATN)
C.
C.    ******************************************************************
C.    *                                                                *
C     *      ROUTINE MULTIPLIES ROTATION MATRIX RMAT BY THE            *
C     *       IROT TH MATRIX IN THE GEANT ROTATION MATRIX BANK         *
C     *     TO PRODUCE THE MATRIX RMATN SUCH THAT RMATN IS             *
C     *       THE RESULT OF FIRST RMAT THEN IROT.                      *
C     *                                                                *
C.    *    ==>Called by : GTRMUL                                       *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION RMAT(10),RMATN(10)
C.
C.    ------------------------------------------------------------------
C.
      IF(IROT.EQ.0) THEN
         DO 10 I=1,10
            RMATN(I)=RMAT(I)
   10    CONTINUE
      ELSE
         JR=LQ(JROTM-IROT)
         RMATN(10)=1.
         IF(RMAT(10).EQ.0.0) THEN
            DO 20 I=1,9
               RMATN(I)=Q(JR+I)
   20       CONTINUE
         ELSE
C
            DO 30 KO=1,9
               RMATN(KO)=0.0
C
               KO1=(KO-1)/3+1
               KO2=KO-(KO1-1)*3
C
               DO 30 IN=1,3
                  KI1=KO2+(IN-1)*3
                  KI2=IN+(KO1-1)*3
                  RMATN(KO)=RMATN(KO)+RMAT(KI1)*Q(JR+KI2)
   30       CONTINUE
         ENDIF
      ENDIF
      END
