*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:16  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GPPARX (NUMB  )
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print particle definition JPART               *
C.    *       NUMB     Particle number                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  S.Giani *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      CHARACTER*32 CHLINE,CHSTRI(50)
      CHARACTER*20 NAME
      DIMENSION ITITLE(5)
C.
C.    ------------------------------------------------------------------
C.
      IF (JPART.LE.0) GO TO 999
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,10000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = IQ(JPART-2)
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
      DO 20 I=N1,N2
         JP = LQ(JPART-I)
         IF (JP.NE.0)THEN
            IOPT = Q(JP+6)
            NL = IQ(JP-1)
            CALL UCOPY(IQ(JP+1),ITITLE,5)
*            NW=MIN(NL,13)
*            NW=NL
            CHLINE='Particle Number='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(I,CHLINE(ILEN:))
            CHSTRI(1)=CHLINE
            CALL UHTOC(IQ(JP+1),4,NAME,20)
            CHLINE='Name='
            ILEN=LENOCC(CHLINE)+1
            CHLINE(ILEN:)=NAME
            CHSTRI(2)=CHLINE
            CHLINE='Iopt='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(IOPT,CHLINE(ILEN:))
            CHSTRI(3)=CHLINE
            CHLINE='Mass='
            ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JP+7),CHLINE(ILEN:))
            WRITE(CHLINE(ILEN:),10200)Q(JP+7)
            CHSTRI(4)=CHLINE
            CHLINE='Charge='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JP+8),CHLINE(ILEN:))
            CHSTRI(5)=CHLINE
            CHLINE='LifeTim='
            ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JP+9),CHLINE(ILEN:))
            WRITE(CHLINE(ILEN:),10300)Q(JP+9)
            CHSTRI(6)=CHLINE
            DO 10 JJ=7,NL
               CHLINE='User='
               ILEN=LENOCC(CHLINE)+1
*       CALL IZRTOC(Q(JP+JJ),CHLINE(ILEN:))
               WRITE(CHLINE(ILEN:),10300)Q(JP+JJ)
               CHSTRI(JJ)=CHLINE
   10       CONTINUE
            CALL IGMESS(NL,CHSTRI,'PRINT','P')
*            WRITE (CHMAIL,1002) I,ITITLE,IOPT,
*     +                          (Q(JP + J),J = 7,NW)
*    5       CALL GMAIL(0,0)
*            IF(NL-NW.GT.0) THEN
*               NS=NW+1
*               NW=MIN(NL,NW+5)
*               WRITE(CHMAIL,1003) (Q(JP + J),J = NS,NW)
*               GO TO 5
*            END IF
         ENDIF
   20 CONTINUE
C
10000 FORMAT ('0',51('='),3X,'Particle Types',3X,50('='))
10100 FORMAT ('0','Part',25X,'Options',8X,'Mass',4X,'Charge'
     +,'    Life time                  User words')
* 1002 FORMAT (' ',I4,1X,5A4,I8,6X,E11.4,F7.0,3X,5(E12.5,2X))
10200 FORMAT (E11.4)
10300 FORMAT (E12.5)
  999 END
