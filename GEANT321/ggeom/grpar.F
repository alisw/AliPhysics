*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:55  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GRPAR(SHAP,NP,P,NDIV,STEP,C0,IN,GP,NSP,SP)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       ROUTINE TO COMPUTE RMIN AND RMAX FOR THE IN TH DIVISION  *
C.    *       OF NDIV DIVISIONS STARTING AT C0 WITH STEP OF STEP       *
C.    *       (FOR ISHAPE 5 & 6) AND TO INSERT THEM INTO THE SHAPE     *
C.    *       PARAMETER LIST FROM P OUTPUTING THE RESULT AS THE        *
C.    *       LIST SP.                                                 *
C.    *                                                                *
C.    *    ==>Called by : not used                                     *
C.    *         Author  A.McPherson  *********                         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcunit.inc"
      DIMENSION P(3),SP(5),GP(5)
C.
C.    ------------------------------------------------------------------
C.
      ISH=SHAP
      IF(ISH.NE.5)GO TO 5
      NSP=3
      SP(3)=P(3)
      SP(2)=C0+STEP*IN
      SP(1)=SP(2)-STEP
      GO TO 99
C
   5  IF(ISH.LT.5.OR.ISH.GT.12) GO TO 90
      IF(ISH.EQ.10) GO TO 90
C
      NSP=NP
      DO 10 I=1,NP
  10  SP(I)=P(I)
C
      IF(ISH.GT.6.AND.ISH.NE.9) GO TO 20
      SP(2)=C0+STEP*IN
      SP(1)=SP(2)-STEP
      GO TO 99
C
   20 CONTINUE
C
      IF(ISH.GT.9) GO TO 30
C
      DR=(GP(3)-GP(2))/NDIV
      SP(3)=GP(2)+DR*IN
      SP(2)=SP(3)-DR
      DR=(GP(5)-GP(4))/NDIV
      SP(5)=GP(4)+DR*IN
      SP(4)=SP(5)-DR
      GO TO 99
C
   30 CONTINUE
C
C               POLYGONS AND POLYCONES.
C
      IPNZ=3
      IF(ISH.EQ.11) IPNZ=4
      NZ=GP(IPNZ)
C
      DO 40 IZ=1,NZ
C
      IPZ=IPNZ+IZ*3-2
      OR=GP(IPZ+1)
      ST=(GP(IPZ+2)-OR)/NDIV
      SP(IPZ+2)=OR+ST*IN
      SP(IPZ+1)=SP(IPZ+2)-ST
C
   40 CONTINUE
C
      GO TO 99
C
  90  CONTINUE
      WRITE(CHMAIL,1000) ISH
      CALL GMAIL(0,0)
C
 1000 FORMAT(' SHAPE #',I5,' AS A DIVISION IN R NOT YET IMPLEMENTED')
  99  RETURN
      END
