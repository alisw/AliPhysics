*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GTRMUL(DX1,RMAT1,DX2,IROT,DXNEW,RMATN)
C.
C.    ******************************************************************
C     *                                                                *
C     *       ROUTINE TO FORM THE TRANSFORMATION RESULTING FROM        *
C     *       THE APPLICATION OF FIRST THE TRANSFORMATION GIVEN        *
C     *       BY THE POINT DX1 AND THE MATRIX RMAT1 FOLLOWED BY        *
C     *       THE APPLICATION OF THE TRANSFORMATION GIVEN BY THE       *
C     *       POINT DX2 AND THE MATRIX IN THE GEANT ROTATION           *
C     *        MATRIX BANKS WITH INDEX NUMBER IROT. THE RESULTING      *
C     *       TRANSFORMATION IS STORED AS THE POINT DXNEW AND          *
C     *       MATRIX RMATN.                                            *
C.    *                                                                *
C.    *    ==>Called by : GFTRAC,GINVOL,GLVOLU,GMEDIV,GMEPOS,GSTRAC,   *
C.    *                   GTMEDI                                       *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
      DIMENSION DX1(3),RMAT1(10),DX2(3),DXNEW(3),RMATN(10)
C.
C.    ------------------------------------------------------------------
C.
      IF(RMAT1(10).EQ.0.0)THEN
         DXNEW(1) = DX1(1) + DX2(1)
         DXNEW(2) = DX1(2) + DX2(2)
         DXNEW(3) = DX1(3) + DX2(3)
         IF(IROT.EQ.0)THEN
            DO 20 I = 1, 10, 2
               RMATN(I) = RMAT1(I)
               RMATN(I+1) = RMAT1(I+1)
   20       CONTINUE
         ELSE
            JR=LQ(JROTM-IROT)
            DO 30 I = 1, 9
               RMATN(I) = Q(JR+I)
   30       CONTINUE
            RMATN(10)=1.
         ENDIF
      ELSE
C
         DXTEM1 = DX2(1)*RMAT1(1) + DX2(2)*RMAT1(4) + DX2(3)*RMAT1(7)
         DXTEM2 = DX2(1)*RMAT1(2) + DX2(2)*RMAT1(5) + DX2(3)*RMAT1(8)
         DXTEM3 = DX2(1)*RMAT1(3) + DX2(2)*RMAT1(6) + DX2(3)*RMAT1(9)
C
         IF(IROT.EQ.0)THEN
            DO 50 I = 1, 10, 2
               RMATN(I) = RMAT1(I)
               RMATN(I+1) = RMAT1(I+1)
   50       CONTINUE
         ELSE
            CALL GRMUL(RMAT1,IROT,RMATN)
         ENDIF
         DXNEW(1) = DX1(1) + DXTEM1
         DXNEW(2) = DX1(2) + DXTEM2
         DXNEW(3) = DX1(3) + DXTEM3
      ENDIF
      END
