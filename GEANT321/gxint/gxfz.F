*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXFZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      FZEBRA control commands                                   *
C.    *                                                                *
C.    *       Authors:   F.Carminati **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcxlun.inc"
#include "geant321/gcflag.inc"
      COMMON/QUEST/IQUEST(100)
      CHARACTER*4  CHOPT
      CHARACTER*80 CHPATH
      CHARACTER*32 CHPATL
      CHARACTER*4  KEYSU
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF (CHPATL.EQ.'FZOPEN') THEN
         CALL KUGETI(LUN)
         CALL GXLUNF(LUN,1,IFREE)
         IF(IFREE.NE.0) GOTO 99
         CALL KUGETC(CHPATH,NCHPAT)
         CALL KUGETC(CHOPT,NCH)
         CALL KUGETI(LRECL)
         CALL CUTOL(CHPATH)
         IF(INDEX(CHOPT,'A').NE.0) THEN
            OPEN(UNIT=LUN,FILE=CHPATH(1:NCHPAT),FORM='FORMATTED',
     +      STATUS='UNKNOWN')
         ELSE
            OPEN(UNIT=LUN,FILE=CHPATH(1:NCHPAT),FORM='UNFORMATTED',
     +      STATUS='UNKNOWN')
         ENDIF
         CALL GOPEN(LUN, CHOPT, LRECL, IER)
         IF(INDEX(CHOPT,'I').NE.0) THEN
            LUNIT(LUN)=1
         ELSE
            LUNIT(LUN)=2
         ENDIF
         IF(NPAR.GT.4) THEN
            CALL KUGETC(CHOPT,NCH)
            IF(INDEX(CHOPT,'I').NE.0)THEN
               LUNIT(LUN)=0
               IDENT = 0
               CALL GFIN(LUN, 'INIT', 1, IDENT, ' ',IER)
               IF(IER.LT.0) THEN
                  WRITE(CHMAIL,12300) IER
12300 FORMAT(' *** GXFZ : Error in GFIN, IER = ',I4)
                  CALL GMAIL(0,0)
               ENDIF
               CALL GCLOSE(LUN,IER)
               LUNIT(LUN)=0
            ELSEIF(INDEX(CHOPT,'O').NE.0)THEN
               LUNIT(LUN)=0
               IDENT=1
               CALL GFOUT(LUN, 'INIT', 1, IDENT, ' ',IER)
               IF(IER.LT.0) THEN
                  WRITE(CHMAIL,12301) IER
12301 FORMAT(' *** GXFZ : Error in GFOUT, IER = ',I4)
                  CALL GMAIL(0,0)
               ENDIF
               CALL GCLOSE(LUN,IER)
               LUNIT(LUN)=0
            ENDIF
         ENDIF
*
      ELSEIF (CHPATL.EQ.'FZIN') THEN
         CALL KUGETI(LUN)
         CALL KUGETC(KEYSU,NCH)
         CALL KUGETI(IDENT)
         CALL GFIN(LUN, KEYSU, 1, IDENT, ' ',IER)
         IF(IER.LT.0) THEN
            WRITE(CHMAIL,12300) IER
            CALL GMAIL(0,0)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'FZOUT') THEN
         CALL KUGETI(LUN)
         CALL KUGETC(KEYSU,NCH)
         CALL KUGETI(IDENT)
         CALL GFOUT(LUN, KEYSU, 1, IDENT, ' ',IER)
         IF(IER.LT.0) THEN
            WRITE(CHMAIL,12301) IER
            CALL GMAIL(0,0)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'FZCLOSE') THEN
         CALL KUGETI(LUN)
         CALL GCLOSE(LUN, IER)
         LUNIT(LUN)=0
*
      ENDIF
*
   99 END
