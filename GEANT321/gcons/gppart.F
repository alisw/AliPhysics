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
      SUBROUTINE GPPART (NUMB  )
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print particle definition JPART               *
C.    *       NUMB     Particle number                                 *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
      DIMENSION ITITLE(5)
C.
C.    ------------------------------------------------------------------
C.
      IF (JPART.LE.0)  GO TO 99
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,1000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = IQ(JPART-2)
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
C
      IF(NUMB.GE.0) THEN
         WRITE (CHMAIL,1001)
         CALL GMAIL(0,1)
      ENDIF
C
      DO 10 I=N1,N2
         JP = LQ(JPART-I)
         IF (JP.NE.0)THEN
            IOPT = Q(JP+6)
            NL = IQ(JP-1)
            CALL UCOPY(IQ(JP+1),ITITLE,5)
            NW=MIN(NL,13)
            WRITE (CHMAIL,1002) I,ITITLE,IOPT,
     +                          (Q(JP + J),J = 7,NW)
    5       CALL GMAIL(0,0)
            IF(NL-NW.GT.0) THEN
               NS=NW+1
               NW=MIN(NL,NW+5)
               WRITE(CHMAIL,1003) (Q(JP + J),J = NS,NW)
               GO TO 5
            END IF
         ENDIF
   10 CONTINUE
C
 1000 FORMAT ('0',51('='),3X,'Particle Types',3X,50('='))
 1001 FORMAT ('0','Part',25X,'Options',8X,'Mass',4X,'Charge'
     +,'    Life time                  User words')
 1002 FORMAT (' ',I4,1X,5A4,I8,6X,E11.4,F7.0,3X,5(E12.5,2X))
 1003 FORMAT (61X,5(E12.5,2X))
  99  RETURN
      END
