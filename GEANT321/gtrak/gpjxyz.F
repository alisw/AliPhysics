*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GPJXYZ(IT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints JXYZ bank for track number IT                     *
C.    *       If IT=0 Prints all track banks                           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
C.
C.    ------------------------------------------------------------------
C.
      I1=1
      I2=NTRACK
      IF(IT.NE.0)I1=IT
      IF(IT.NE.0)I2=IT
      IF(I2.GT.NTRACK)GO TO 99
      IF(I1.LT.1)GO TO 99
      IF(JKINE.LE.0)GO TO 99
      IF(JXYZ.LE.0)GO TO 99
C
      DO 90 I=I1,I2
         JK=LQ(JKINE-I)
         IF(JK.NE.0)THEN
            IPART=Q(JK+5)
            JPA=LQ(JPART-IPART)
            WRITE(CHMAIL,1000)I,(IQ(JPA+L),L=1,4)
            CALL GMAIL(0,0)
            WRITE(CHMAIL,1100)
            CALL GMAIL(0,0)
C
            JX=LQ(JXYZ-I)
            IPOINT=JX+3
            IP=0
  10        NPOINT=Q(IPOINT)
            IF(NPOINT.GT.0)THEN
               IPT=Q(IPOINT+1)
               DO 20 K=1,NPOINT
                  KP=IPOINT+3*K-2
                  IP=IP+1
                  WRITE(CHMAIL,2000)IP,IPT,(Q(KP+L),L=1,3)
                  CALL GMAIL(0,0)
  20           CONTINUE
C
               IPOINT=IPOINT+3.*Q(IPOINT)+2.
               GO TO 10
            ENDIF
         ENDIF
  90  CONTINUE
C
 1000 FORMAT(3X,'JXYZ',3X,'TRACK NUMBER',I5,2X,4A4)
 1100 FORMAT(6X,'POINT',6X,'PART',6X,'X',10X,'Y',10X,'Z')
 2000 FORMAT(4X,I6,6X,I4,2X,F10.4,1X,F10.4,1X,F10.4)
C
  99  RETURN
      END
