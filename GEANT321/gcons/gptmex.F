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
      SUBROUTINE GPTMEX (NUMB  )
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Routine to print tracking media data structure JTMED     *
C.    *       NUMB     medium number                                   *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GPRINT                               *
C.    *       Author    R.Brun S.Giani  *********                      *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcnum.inc"
      CHARACTER*32 CHLINE,CHSTRI(50)
      CHARACTER*20 NAME
      DIMENSION ITITLE(5)
C.
C.    ------------------------------------------------------------------
C.
      IF (JTMED.LE.0) GO TO 999
C
      IF (NUMB  .EQ.0)THEN
         WRITE (CHMAIL,10000)
         CALL GMAIL(0,0)
         N1     = 1
         N2     = NTMED
      ELSE
         N1     = ABS(NUMB)
         N2     = ABS(NUMB)
      ENDIF
*      IF(NUMB.GE.0) THEN
*         WRITE (CHMAIL,1001)
*         CALL GMAIL(0,1)
*      ENDIF
C
      DO 10 I=N1,N2
         JTM = LQ(JTMED-I)
         IF (JTM.NE.0)THEN
            IMAT = Q(JTM+6)
            ISVOL = Q(JTM+7)
            IFIELD = Q(JTM+8)
            CALL UCOPY(IQ(JTM+1),ITITLE,5)
            CHLINE='T. Med. Number='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(I,CHLINE(ILEN:))
            CHSTRI(1)=CHLINE
            CALL UHTOC(IQ(JTM+1),4,NAME,20)
            CHLINE='Name='
            ILEN=LENOCC(CHLINE)+1
            CHLINE(ILEN:)=NAME
            CHSTRI(2)=CHLINE
            CHLINE='Material='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(IMAT,CHLINE(ILEN:))
            CHSTRI(3)=CHLINE
            CHLINE='Isvol='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(ISVOL,CHLINE(ILEN:))
            CHSTRI(4)=CHLINE
            CHLINE='Ifield='
            ILEN=LENOCC(CHLINE)+1
            CALL IZITOC(IFIELD,CHLINE(ILEN:))
            CHSTRI(5)=CHLINE
            CHLINE='Fieldm='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JTM+9),CHLINE(ILEN:))
            CHSTRI(6)=CHLINE
            CHLINE='Tmaxfd='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JTM+10),CHLINE(ILEN:))
            CHSTRI(7)=CHLINE
            CHLINE='Stepmax='
            ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JTM+11),CHLINE(ILEN:))
            WRITE(CHLINE(ILEN:),10200)Q(JTM+11)
            CHSTRI(8)=CHLINE
            CHLINE='Deemax='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JTM+12),CHLINE(ILEN:))
            CHSTRI(9)=CHLINE
            CHLINE='Epsil='
            ILEN=LENOCC(CHLINE)+1
*      CALL IZRTOC(Q(JTM+13),CHLINE(ILEN:))
            WRITE(CHLINE(ILEN:),10200)Q(JTM+13)
            CHSTRI(10)=CHLINE
            CHLINE='Stmin='
            ILEN=LENOCC(CHLINE)+1
            CALL IZRTOC(Q(JTM+14),CHLINE(ILEN:))
            CHSTRI(11)=CHLINE
            CALL IGMESS(11,CHSTRI,'PRINT','P')
 
*            WRITE(CHMAIL,1002)I,ITITLE,IMAT,ISVOL,IFIELD,
*     +                        (Q(JTM+J),J=9,14)
*            CALL GMAIL(0,0)
         ENDIF
   10 CONTINUE
C
10000 FORMAT ('0',51('='),3X,'TRACKING MEDIA',3X,50('='))
10100 FORMAT ('0','TMED',26X,'MATERIAL ','ISVOL',' IFIELD  FIELDM'
     +, '  TMAXFD','  STEMAX','    DEEMAX','   EPSIL','   STMIN')
10200 FORMAT (E12.5)
  999 END
