*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:48  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFCORD(NIN,CLOW,CHIGH,CORD,ITYPE,ICON)
C.
C.    *****************************************************************
C.    *                                                               *
C.    *    ROUTINE TO ORDER THE VALUES IN CLOW AND CHIGH INTO CORD    *
C.    *    WHILST STORING A 1 IN ITYPE IF THE VALUE IS FROM CLOW      *
C.    *    AND A 2 IF FROM CHIGH. THE SEQUENCE NUMBER OF THE VALUE    *
C.    *    IS STORED IN ICON. THERE ARE NIN VALUES IN EACH OF CLOW    *
C.    *    AND CHIGH.                                                 *
C.    *                                                               *
C.    *    ==>Called by : GGORD                                       *
C.    *         Author  A.McPherson  *********                        *
C.    *                                                               *
C.    *****************************************************************
C.
      DIMENSION CLOW(2),CHIGH(2),CORD(2),ITYPE(2),ICON(2)
C.
C.       -----------------------------------------------------
C.
      IC=0
      DO 60 IN=1,NIN
C
      ITRY=1
      CTRY=CLOW(IN)
C
      IPL=1
      IF(IC.EQ.0) GO TO 50
C
   10 CONTINUE
C
      DO 20 I=1,IC
      IF(CORD(I).LE.CTRY) GO TO 20
      IPL=I
      GO TO 30
   20 CONTINUE
C
      IPL=IC+1
      GO TO 50
C
   30 CONTINUE
      DO 40 I=IPL,IC
C
      IC1=IC+IPL-I
      IC2=IC1+1
      CORD(IC2)=CORD(IC1)
      ICON(IC2)=ICON(IC1)
      ITYPE(IC2)=ITYPE(IC1)
C
   40 CONTINUE
C
   50 CONTINUE
C
      IC=IC+1
      CORD(IPL)=CTRY
      ITYPE(IPL)=ITRY
      ICON(IPL)=IN
C
      IF(ITRY.NE.1) GO TO 60
      ITRY=2
      CTRY=CHIGH(IN)
      GO TO 10
C
   60 CONTINUE
C
      END
