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
      SUBROUTINE GPKINX(IT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Prints kinematics bank JKINE for track number IT         *
C.    *       If IT=0 Prints all kinematics banks                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  S.Giani *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcnum.inc"
#include "geant321/gcunit.inc"
      CHARACTER*20 NAME
      CHARACTER*32 CHLINE,CHSTRI(50)
C.
C.    ------------------------------------------------------------------
C.
*      WRITE(CHMAIL,1000)
*      CALL GMAIL(1,0)
C
      I1=1
      I2=NTRACK
      IF(IT.NE.0)I1=IT
      IF(IT.NE.0)I2=IT
      IF(I2.GT.NTRACK)GO TO 999
      IF(I1.LT.1)GO TO 999
      IF(JKINE.LE.0)GO TO 999
C
      DO 40 ITRACK=I1,I2
         JK=LQ(JKINE-ITRACK)
         IF(JK.NE.0)THEN
            IPART=Q(JK+5)
            JPA=LQ(JPART-IPART)
            NL=IQ(JK-2)
            NG=Q(JK+7)+7
            NW=MIN(NG,16)
            CHLINE='Kine='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(ITRACK,CHLINE(ILEN:))
            CHSTRI(1)=CHLINE
            CALL UHTOC(IQ(JPA+1),4,NAME,20)
            CHLINE='Name='
            ILEN=LENOCC(CHLINE)+1
            CHLINE(ILEN:)=NAME
            CHSTRI(2)=CHLINE
            CHLINE='Type='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(IPART,CHLINE(ILEN:))
            CHSTRI(3)=CHLINE
            CHLINE='Px='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JK+1),CHLINE(ILEN:))
            CHSTRI(4)=CHLINE
            CHLINE='Py='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JK+2),CHLINE(ILEN:))
            CHSTRI(5)=CHLINE
            CHLINE='Pz='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JK+3),CHLINE(ILEN:))
            CHSTRI(6)=CHLINE
            CHLINE='Energy='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JK+4),CHLINE(ILEN:))
            CHSTRI(7)=CHLINE
            CHLINE='Vorigin='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JK+6),CHLINE(ILEN:))
*      WRITE(CHLINE(ILEN:),1002)Q(JK+6)
            CHSTRI(8)=CHLINE
            DO 10 JJ=8,NW
               CHLINE='Vgener='
               ILEN=LENOCC(CHLINE)+1
               CALL IZITOC(INT(Q(JK+JJ)),CHLINE(ILEN:))
*       WRITE(CHLINE(ILEN:),1002)INT(Q(JK+JJ)
               CHSTRI(1+JJ)=CHLINE
   10       CONTINUE
            NLINE=8
            IF(NW.GE.8)NLINE=8+(NW-8+1)
            CALL IGMESS(NLINE,CHSTRI,'PICK','P')
*            WRITE(CHMAIL,2000)ITRACK,IQ(JPA+1),IQ(JPA+2),IQ(JPA+3),IPART
*     +       ,(Q(JK+I),I=1,4),Q(JK+6),(INT(Q(JK+L)),L=8,NW)
*   10       CALL GMAIL(0,0)
   20       CONTINUE
            IF(NG-NW.GT.0) THEN
               NS=NW+1
               NW=MIN(NG,NW+9)
               WRITE(CHMAIL,10200) (INT(Q(JK+L)),L=NS,NW)
               CALL GMAIL(0,0)
               GO TO 20
            END IF
            IF(NL.GT.0) THEN
*               WRITE(CHMAIL,10300)
*               CALL GMAIL(0,0)
*               JUK = LQ(JK-1)
*               NL  = IQ(JUK-1)
*               NW=MIN(NL,9)
*               WRITE(CHMAIL,10400) (Q(JUK+L),L=1,NW)
*   30          CALL GMAIL(0,0)
*               IF(NL-NW.GT.0) THEN
*                  NS=NW+1
*                  NW=MIN(NL,NW+9)
*                  WRITE(CHMAIL,10400) (Q(JUK+L),L=NS,NW)
*                  GO TO 30
*               END IF
            END IF
         ENDIF
   40 CONTINUE
C
10000 FORMAT(3X,'Kine',5X,'Name',10X,'Type',5X,'Px',9X,'Py',
     +       9X,'Pz',7X,'Energy',5X,'Vorigin',2X,'Vgener')
10100 FORMAT(1X,I5,6X,3A4,I5,5(1X,F10.4),9(1X,I4))
10200 FORMAT(84X,9(1X,I4))
10300 FORMAT(62X,'User words')
10400 FORMAT(1X,9(2X,E12.5))
10500 FORMAT(E12.5)
C
  999 END
