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
      SUBROUTINE GPTMED (NUMB  )
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print tracking media data structure JTMED     *
C.    *       NUMB     medium number                                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      DIMENSION ITITLE(5)
C.
C.    ------------------------------------------------------------------
C.
      IF (JTMED.LE.0)  GO TO 99
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,1000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = NTMED
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
      IF(NUMB.GE.0) THEN
         WRITE (CHMAIL,1001)
         CALL GMAIL(0,1)
      ENDIF
C
      DO 10 I=N1,N2
         JTM = LQ(JTMED-I)
         IF (JTM.NE.0)THEN
            IMAT = Q(JTM+6)
            ISVOL = Q(JTM+7)
            IFIELD = Q(JTM+8)
            CALL UCOPY(IQ(JTM+1),ITITLE,5)
            WRITE(CHMAIL,1002)I,ITITLE,IMAT,ISVOL,IFIELD,
     +                        (Q(JTM+J),J=9,14)
            CALL GMAIL(0,0)
         ENDIF
   10 CONTINUE
C
 1000 FORMAT ('0',51('='),3X,'TRACKING MEDIA',3X,50('='))
 1001 FORMAT ('0','TMED',26X,'MATERIAL ','ISVOL',' IFIELD  FIELDM'
     +, '  TMAXFD','  STEMAX','    DEEMAX','   EPSIL','   STMIN')
 1002 FORMAT (' ',I6,1X,5A4,I8,I8,I6,4X,F6.2,2X,F6.2,G10.3,3F8.3)
  99  RETURN
      END
