*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.19  by  S.Giani
*-- Author :
      SUBROUTINE GPMATE (NUMB)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print material data structures JMATE          *
C.    *       NUMB     Material number                                 *
C.    *                                                                *
C.    *      Changed by S.Egli at 8.5.90: also show mixture contents   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      CHARACTER CHMIXT*17
C.
C.    ------------------------------------------------------------------
C.
      IF (JMATE.LE.0) GO TO 999
      IF (NUMB .EQ.0) THEN
         WRITE (CHMAIL,10000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = NMATE
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
      IF(NUMB.GE.0) THEN
         WRITE (CHMAIL,10100)
         CALL GMAIL(0,1)
      ENDIF
C
      DO 20 I=N1,N2
         JMA = LQ(JMATE-I)
         IF (JMA.LE.0) GO TO 20
C
         NMIXT=Q(JMA+11)
         CHMIXT=' '
         IF(NMIXT.GT.1) CHMIXT='   A      Z     W'
         WRITE (CHMAIL,10200) I,(Q(JMA + J),J = 1,10),NMIXT,CHMIXT
         CALL GMAIL(0,0)
         IF(NMIXT.GT.1)THEN
            JMX=LQ(JMA-5)
            DO 10 J=1,NMIXT
               WRITE(CHMAIL,10300)Q(JMX+J),Q(JMX+NMIXT+J),
     +         Q(JMX+2*NMIXT+J)
               CALL GMAIL(0,0)
   10       CONTINUE
         ENDIF
   20 CONTINUE
C
10000 FORMAT ('0',51('='),5X,'MATERIALS',6X,50('='))
10100 FORMAT ('0','MATERIAL',27X,'A',9X,'Z',5X,'DENSITY'
     +,2X,'RADIAT L',2X,'ABSORP L',' NMIXT')
10200 FORMAT (' ',I8,1X,5A4,3F10.3,2E10.3,I4,2X,A17)
10300 FORMAT (' ',85X,2F7.2,F7.3)
  999 END
