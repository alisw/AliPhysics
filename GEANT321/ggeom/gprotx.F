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
      SUBROUTINE GPROTX(NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print rotation matrices structure JROTM       *
C.    *       NUMB     Rotation matrix number                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *         Author  R.Brun  S.Giani *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*32 CHLINE,CHSTRI(50)
C.
C.    ------------------------------------------------------------------
C.
      IF (JROTM.LE.0) GO TO 999
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,10000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = IQ(JROTM-2)
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
C
*      IF(NUMB.GE.0) THEN
*         WRITE (CHMAIL,1001)
*         CALL GMAIL(0,1)
*      ENDIF
C
      DO 10 I=N1,N2
         JR = LQ(JROTM-I)
         IF (JR.LE.0) GO TO 10
C
         IFLAG = Q(JR+10)
         CHLINE='Matrix Number='
         ILEN=LENOCC(CHLINE)+1
         CALL IZITOC(I,CHLINE(ILEN:))
         CHSTRI(1)=CHLINE
         CHLINE='RM(1,1)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+1),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+1)
         CHSTRI(2)=CHLINE
         CHLINE='RM(2,1)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+2),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+2)
         CHSTRI(3)=CHLINE
         CHLINE='RM(3,1)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+3),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+3)
         CHSTRI(4)=CHLINE
         CHLINE='RM(1,2)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+4),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+4)
         CHSTRI(5)=CHLINE
         CHLINE='RM(2,2)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+5),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+5)
         CHSTRI(6)=CHLINE
         CHLINE='RM(3,2)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+6),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+6)
         CHSTRI(7)=CHLINE
         CHLINE='RM(1,3)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+7),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+7)
         CHSTRI(8)=CHLINE
         CHLINE='RM(2,3)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+8),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+8)
         CHSTRI(9)=CHLINE
         CHLINE='RM(3,3)='
         ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JR+9),CHLINE(ILEN:))
         WRITE(CHLINE(ILEN:),10200)Q(JR+9)
         CHSTRI(10)=CHLINE
         CHLINE='Iflag='
         ILEN=LENOCC(CHLINE)+1
         CALL IZITOC(IFLAG,CHLINE(ILEN:))
         CHSTRI(11)=CHLINE
         CHLINE='Theta1='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+11),CHLINE(ILEN:))
         CHSTRI(12)=CHLINE
         CHLINE='Phi1='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+12),CHLINE(ILEN:))
         CHSTRI(13)=CHLINE
         CHLINE='Theta2='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+13),CHLINE(ILEN:))
         CHSTRI(14)=CHLINE
         CHLINE='Phi2='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+14),CHLINE(ILEN:))
         CHSTRI(15)=CHLINE
         CHLINE='Theta3='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+15),CHLINE(ILEN:))
         CHSTRI(16)=CHLINE
         CHLINE='Phi3='
         ILEN=LENOCC(CHLINE)+1
         CALL IZRTOC(Q(JR+16),CHLINE(ILEN:))
         CHSTRI(17)=CHLINE
         CALL IGMESS(17,CHSTRI,'PRINT','P')
*      WRITE (CHMAIL,1002) I,(Q(JR+J),J = 1,9),IFLAG,(Q(JR+J),J = 11,16)
*      CALL GMAIL(0,0)
   10 CONTINUE
C
10000 FORMAT ('0',51('='),1X,'ROTATION MATRICES',2X,50('='))
10100 FORMAT ('0','MATRIX',1X,'<',17('-'),' ROTATION MATRIX',1X,17('-')
     +,'>',' FLAG','  THET1','   PHI1','  THET2','   PHI2','  THET3'
     +,'   PHI3')
* 1002 FORMAT (' ',I6,9F6.3,I5,6F7.2)
10200 FORMAT (E12.5)
  999 END
