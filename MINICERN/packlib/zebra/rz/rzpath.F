*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:26  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZPATH(CHPATH)
************************************************************************
*
*       RZPATH will decode the character string CHPATH into a character
*        array CHPAT(1),,CHPAT(NLPAT)
*        where NLPAT is the number of levels starting from the root directory
*
* Input:
*   *CHPATH* Character variable  specifying the pathname
*
* Called by RZCDIR,RZCOPY,RZLDIR,RZNDIR
*
*  Author  : R.Brun DD/US/PD
*  Written : 01.04.86
*  Last mod: 15.11.91 J.Shiers Support Unix-style .. as well as \
*
************************************************************************
#include "zebra/zbcdch.inc"
                     CHARACTER*1  BSLASH,KTILDE
#include "zebra/zunit.inc"
#include "zebra/rzdir.inc"
      CHARACTER*(*) CHPATH
      CHARACTER*1 CH1
      CHARACTER*2 CH2
*
*-----------------------------------------------------------------------
*
      BSLASH=CQALLC(61:61)
      KTILDE=CQALLC(94:94)
      NCHP=LEN(CHPATH)
      NLPAT=0
  10  IF(CHPATH(NCHP:NCHP).EQ.' ')THEN
         NCHP=NCHP-1
         IF(NCHP.GT.0)GO TO 10
         NLPAT=NLCDIR
         DO 20 I=1,NLCDIR
            CHPAT(I)=CHCDIR(I)
  20     CONTINUE
         GO TO 99
      ENDIF
*
      IS1=1
  30  IF(CHPATH(IS1:IS1).EQ.' ')THEN
         IS1=IS1+1
         GO TO 30
      ENDIF
*
      CH1=CHPATH(IS1:IS1)
      CH2=CHPATH(IS1:IS1+1)
      IF(CH1.EQ.'/')THEN
         IF(CHPATH(IS1+1:IS1+1).EQ.'/')THEN
*
*           Case where Path starts with //Topdirectory
*
            IS=IS1+2
            IF(IS.GT.NCHP)GO TO 99
  40        IF(CHPATH(IS:IS).EQ.'/')THEN
               IF(IS.EQ.IS1+2)GO TO 90
               NLPAT=1
               CHPAT(1)=CHPATH(IS1+2:IS-1)
               IS1=IS+1
               IS=IS1
               GO TO 50
            ELSE
               IS=IS+1
               IF(IS.LT.NCHP)GO TO 40
               NLPAT=1
               CHPAT(1)=CHPATH(IS1+2:IS)
               GO TO 99
            ENDIF
         ENDIF
*
*           Case where Path starts just below topdirectory /a/b/c
*
         IF(CHPATH(IS1+1:IS1+1).EQ.BSLASH)GO TO 90
         IF(CHPATH(IS1+1:IS1+1).EQ.KTILDE)GO TO 90
         NLPAT=1
         CHPAT(1)=CHCDIR(1)
         IS=IS1+1
         IS1=IS
*
  50     IF(IS.EQ.NCHP)THEN
            IF(CHPATH(IS1:IS).NE.'..'.AND.
     +         CHPATH(IS1:IS).NE.BSLASH) THEN
               NLPAT=NLPAT+1
               IF(NLPAT.GT.NLPATM)GO TO 90
               CHPAT(NLPAT)=CHPATH(IS1:IS)
            ELSE
               NLPAT = NLPAT -1
            ENDIF
            GO TO 99
         ELSE
            IF(CHPATH(IS:IS).EQ.'/')THEN
               IF(NLPAT.GT.NLPATM)GO TO 90
               IF(CHPATH(IS1:IS-1).NE.'..'.AND.
     +            CHPATH(IS1:IS-1).NE.BSLASH) THEN
                  NLPAT=NLPAT+1
                  CHPAT(NLPAT)=CHPATH(IS1:IS-1)
               ELSE
                  NLPAT = NLPAT - 1
               ENDIF
               IS1=IS+1
            ENDIF
            IS=IS+1
            GO TO 50
         ENDIF
      ENDIF
*
*          Case where Path is relative to the naming directory
*
      IF(CH1.EQ.KTILDE)THEN
         NLPAT=NLNDIR
         DO 60 I=1,NLNDIR
            CHPAT(I)=CHNDIR(I)
  60     CONTINUE
         IF(IS1.EQ.NCHP)GO TO 99
         IS1=IS1+1
         CH1=CHPATH(IS1:IS1)
         GO TO 75
      ENDIF
*
*          Case where Path is relative to the working directory
*
      DO 70 I=1,NLCDIR
         CHPAT(I)=CHCDIR(I)
  70  CONTINUE
      NLPAT=NLCDIR
*
  75  IF(CH1.EQ.BSLASH)THEN
         NLPAT=NLPAT-1
         IF(NLPAT.EQ.0)NLPAT=1
         IF(IS1.EQ.NCHP)GO TO 99
         IS1=IS1+1
         CH1=CHPATH(IS1:IS1)
         GO TO 75
      ENDIF
*
      IS=IS1
  76  IF(CH2.EQ.'..')THEN
         NLPAT=NLPAT-1
         IF(NLPAT.EQ.0)NLPAT=1
         IF(IS1+1.EQ.NCHP)GO TO 99
         IF(CHPATH(IS1+2:IS1+2).NE.'/') GOTO 90
         IS =IS1
         IS1=IS1+3
         CH2=CHPATH(IS1:IS1+1)
         GO TO 76
      ENDIF
*
  80  IF(IS.EQ.NCHP)THEN
         NLPAT=NLPAT+1
         IF(NLPAT.GT.NLPATM)GO TO 90
         CHPAT(NLPAT)=CHPATH(IS1:IS)
         GO TO 99
      ELSE
         IF(CHPATH(IS:IS).EQ.'/')THEN
            IF(IS.GT.IS1)THEN
               NLPAT=NLPAT+1
               IF(NLPAT.GT.NLPATM)GO TO 90
               CHPAT(NLPAT)=CHPATH(IS1:IS-1)
            ENDIF
            IS1=IS+1
         ENDIF
         IS=IS+1
         GO TO 80
      ENDIF
*
*           Error in Pathname
*
  90  NLPAT=0
*
  99  RETURN
      END
