*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPAFF(CH,NL,CHPATH)
************************************************************************
*
*       RZPAFF reconstructs into CHPATH the path-name stored into the
*       character array CH.
*       where NL=number of levels starting from the root directory
*
* Input:
*    CH  Character array of NL elements
* Output:
*   *CHPATH* Character variable  specifying the pathname
*
* Called by RZCDIR,RZFDIR,RZLDIR,RZLLOK,RZNDIR
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 24.06.92 JDS - IQUEST(1) /= 0 in case of problems
*
************************************************************************
      CHARACTER*(*) CHPATH,CH(*)
*     CHARACTER*128 CHTEMP
      CHARACTER*255 CHTEMP
      CHARACTER*16  CHL
#include "zebra/quest.inc"
*
*-----------------------------------------------------------------------
*
      MAXLEN=LEN(CHPATH)
*     IF(MAXLEN.GT.110)MAXLEN=110
      IF(MAXLEN.GT.255)MAXLEN=255
      IQUEST(1) = 0
*
      CHPATH='//'//CH(1)
      LENG=LENOCC(CHPATH)
      IF(LENG.EQ.2) THEN
        CHPATH='//HOME'
        LENG=6
      ENDIF
      IF(NL.EQ.1) GOTO 99
*
      DO 20 I=2,NL
         CHL=CH(I)
         NMAX=LENOCC(CHL)
         IF(NMAX.EQ.0) THEN
            IQUEST(1) = 1
            GOTO 99
         ENDIF
         IF(LENG+NMAX.GT.MAXLEN)NMAX=MAXLEN-LENG
         CHTEMP=CHPATH(1:LENG)//'/'//CHL(1:NMAX)
         CHPATH=CHTEMP
         LENG=LENG+NMAX+1
         IF(LENG.EQ.MAXLEN) THEN
            IQUEST(1) = 2
            GOTO 99
         ENDIF
  20  CONTINUE
*
  99  RETURN
      END
