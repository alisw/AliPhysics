*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:18  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GPKINE(IT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints kinematics bank JKINE for track number IT         *
C.    *       If IT=0 Prints all kinematics banks                      *
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
      WRITE(CHMAIL,1000)
      CALL GMAIL(1,0)
C
      I1=1
      I2=NTRACK
      IF(IT.NE.0)I1=IT
      IF(IT.NE.0)I2=IT
      IF(I2.GT.NTRACK)GO TO 99
      IF(I1.LT.1)GO TO 99
      IF(JKINE.LE.0)GO TO 99
C
      DO 20 ITRACK=I1,I2
         JK=LQ(JKINE-ITRACK)
         IF(JK.NE.0)THEN
            IPART=Q(JK+5)
            JPA=LQ(JPART-IPART)
            NL=IQ(JK-2)
            NG=Q(JK+7)+7
            NW=MIN(NG,16)
            WRITE(CHMAIL,2000)ITRACK,IQ(JPA+1),IQ(JPA+2),IQ(JPA+3),IPART
     +       ,(Q(JK+I),I=1,4),Q(JK+6),(INT(Q(JK+L)),L=8,NW)
   10       CALL GMAIL(0,0)
            IF(NG-NW.GT.0) THEN
               NS=NW+1
               NW=MIN(NG,NW+9)
               WRITE(CHMAIL,3000) (INT(Q(JK+L)),L=NS,NW)
               GO TO 10
            END IF
            IF(NL.GT.0) THEN
               JUK = LQ(JK-1)
               IF (JUK.GT.0) THEN
                  WRITE(CHMAIL,4000)
                  CALL GMAIL(0,0)
                  ND  = IQ(JUK-1)
                  NW=MIN(ND,9)
                  WRITE(CHMAIL,5000) (Q(JUK+L),L=1,NW)
   15             CALL GMAIL(0,0)
                  IF(ND-NW.GT.0) THEN
                     NS=NW+1
                     NW=MIN(ND,NW+9)
                     WRITE(CHMAIL,5000) (Q(JUK+L),L=NS,NW)
                     GO TO 15
                  END IF
               END IF
            END IF
         ENDIF
   20 CONTINUE
C
 1000 FORMAT(3X,'Kine',5X,'Name',10X,'Type',5X,'Px',9X,'Py',
     +       9X,'Pz',7X,'Energy',5X,'Vorigin',2X,'Vgener')
 2000 FORMAT(1X,I5,6X,3A4,I5,5(1X,F10.4),9(1X,I4))
 3000 FORMAT(84X,9(1X,I4))
 4000 FORMAT(62X,'User words')
 5000 FORMAT(1X,9(2X,E12.5))
C
  99  RETURN
      END
