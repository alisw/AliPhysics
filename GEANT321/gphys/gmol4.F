*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:26  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GMOL4(Y,X,VAL,ARG,EPS,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    * ROUTINE TAKEN FROM IBM SCIENTIFIC SUBROUTINE PACKAGE           *
C.    *                                                                *
C.    * 4 POINT CONTINUED FRACTION INTERPOLATION                       *
C.    * Y=INTERPOLATED VALUE FOR THE ARGUMENT X .                      *
C.    * VAL=VALUE ARRAY .                                              *
C.    * ARG=ARGUMENT ARRAY .                                           *
C.    * EPS=DESIRED ACCURACY .                                         *
C.    * OUTPUT ERROR PARAMETER IER = 0 ACCURACY O.K.                   *
C.    *                            = 1 ACCURACY CAN NOT BE TESTED      *
C.    *                                IN 4TH ORDER INTERPOLATION .    *
C.    *                            = 2 TWO IDENTICAL ELEMENTS IN THE   *
C.    *                                ARGUMENT ARRAY .                *
C.    *                                                                *
C.    *    ==>Called by : GMOLIE                                       *
C.    *         Author M.S. Dixit NRCC Ottawa    *********             *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION ARG(4),VAL(4)
C.
C.    ------------------------------------------------------------------
C.
      IER=1
      Y=VAL(1)
      P2=1.
      P3=Y
      Q2=0.
      Q3=1.
      DO 16 I=2,4
         II=0
         P1=P2
         P2=P3
         Q1=Q2
         Q2=Q3
         Z=Y
         JEND=I-1
   3     AUX=VAL(I)
         DO 10 J=1,JEND
            H=VAL(I)-VAL(J)
            IF(ABS(H).GT.1.E-6*ABS(VAL(I)))GO TO 9
            IF(ARG(I).EQ.ARG(J))GO TO 17
            IF(J.LT.JEND)GO TO 8
            II=II+1
            III=I+II
            IF(III.GT.4)GO TO 19
            VAL(I)=VAL(III)
            VAL(III)=AUX
            AUX=ARG(I)
            ARG(I)=ARG(III)
            ARG(III)=AUX
            GO TO 3
   8        VAL(I)=1.E36
            GO TO 10
   9        VAL(I)=(ARG(I)-ARG(J))/H
  10     CONTINUE
         P3=VAL(I)*P2+(X-ARG(I-1))*P1
         Q3=VAL(I)*Q2+(X-ARG(I-1))*Q1
         IF(Q3.NE.0.)THEN
            Y=P3/Q3
         ELSE
            Y=1.E36
         ENDIF
         DELT=ABS(Z-Y)
         IF(DELT.LE.EPS)GO TO 19
  16  CONTINUE
      RETURN
  17  IER=2
      RETURN
  19  IER=0
      END
