*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXRZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      RZEBRA control commands                                   *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                  P.Zanarini  **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcxlun.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcrz.inc"
      COMMON/QUEST/IQUEST(100)
      CHARACTER*20 CHDIR,CHOPT
      CHARACTER*64 CHPATH
      CHARACTER*64 CHKEY
      CHARACTER*8 CHTAG(20)
      DIMENSION KEYRZ(4)
      CHARACTER*32 CHPATL
      CHARACTER*(*) BSLASH
#if defined(CERNLIB_BSLASH)
      PARAMETER (BSLASH='\\')
#endif
#if !defined(CERNLIB_BSLASH)
      PARAMETER (BSLASH='\')
#endif
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF (CHPATL.EQ.'PQUEST') THEN
         CALL KUGETI(IQ1)
         CALL KUGETI(IQ2)
         IF (NPAR.EQ.1) IQ2=IQ1
         DO 10 I=IQ1,IQ2
            WRITE(CHMAIL,10000)I,IQUEST(I),IQUEST(I)
            CALL GMAIL(0,0)
   10    CONTINUE
10000 FORMAT(' IQUEST(',I2,')=',I10,2X,Z8)
*
      ELSEIF (CHPATL.EQ.'FILE') THEN
         CALL KUGETI(LUN)
         CALL GXLUNF(LUN,1,IFREE)
         IF(IFREE.NE.0)GO TO 999
         CALL KUGETC(CHPATH,NCH)
         CALL KUGETC(CHOPT,NCH)
         CALL GRFILE(LUN,CHPATH,CHOPT)
         IF(INDEX(CHOPT,'I').NE.0.OR.INDEX(CHOPT,'O').NE.0)THEN
            LUNIT(LUN)=0
         ELSE
            LUNIT(LUN)=5
         ENDIF
*
      ELSEIF (CHPATL.EQ.'REND') THEN
         CALL KUGETI(IRLUN)
         CALL GREND(IRLUN)
*
      ELSEIF (CHPATL.EQ.'RZSAVE') THEN
         CALL RZSAVE
*
      ELSEIF (CHPATL.EQ.'MDIR') THEN
         CALL HCDIR(CHDIR,'R')
         IF(CHDIR(3:6).EQ.'PAWC')THEN
            IPAW=1
         ELSE
            IPAW=0
         ENDIF
         CALL KUGETC(CHDIR,NCH)
         DO 20 I=1,NCH
            IF(CHDIR(I:I).EQ.'/'.OR.CHDIR(I:I).EQ.BSLASH)THEN
               WRITE(CHMAIL,10100)
10100       FORMAT('Directory name cannot contain slash')
               CALL GMAIL(0,0)
               GOTO 999
            ENDIF
   20    CONTINUE
         IF(IPAW.NE.0)THEN
            CALL HMDIR(CHDIR,' ')
         ELSE
            CALL KUGETI(NWKEY)
            CALL KUGETC(CHOPT,NCH)
            DO 30 I=1,NWKEY
               CALL KUGETS(CHTAG(I),NCH)
   30       CONTINUE
            CALL RZMDIR(CHDIR,NWKEY,CHOPT,CHTAG)
         ENDIF
*
      ELSEIF (CHPATL.EQ.'CDIR') THEN
         IF(NPAR.EQ.0)THEN
            CHOPT='P'
            CHPATH=' '
         ELSE
            CALL KUGETC(CHPATH,NCH)
            CHOPT=' '
            CALL KUGETC(CHOPT,NCH)
         ENDIF
         CALL HCDIR(CHPATH,CHOPT)
*
      ELSEIF (CHPATL.EQ.'IN') THEN
         CALL KUGETC(CHKEY,NCH)
         CALL KUGETI(IDVERS)
         CALL KUGETC(CHOPT,NCH)
         CALL GRIN(CHKEY,IDVERS,CHOPT)
*
      ELSEIF (CHPATL.EQ.'OUT') THEN
         CALL KUGETC(CHKEY,NCH)
         CALL KUGETI(IDVERS)
         CALL KUGETC(CHOPT,NCH)
         CALL GROUT(CHKEY,IDVERS,CHOPT)
*
      ELSEIF (CHPATL.EQ.'LDIR') THEN
         CALL KUGETC(CHPATH,NCH)
         CHOPT='R'
         CALL KUGETC(CHOPT,NCH)
         IF(NCH.GT.0)CHOPT=CHOPT(1:NCH)//'R'
         CALL HLDIR(CHPATH,CHOPT)
*
      ELSEIF (CHPATL.EQ.'PURGE') THEN
         CALL KUGETI(NKEEP)
         CALL RZPURG(NKEEP)
*
      ELSEIF (CHPATL.EQ.'SCR') THEN
         CALL KUGETC(CHKEY,NCH)
         CALL KUGETI(IDVERS)
         CALL UCTOH(CHKEY,KEYRZ(1),4,4)
         KEYRZ(2)=IDVERS
         ICYCLE=9999
         CHOPT=' '
         CALL RZDELK(KEYRZ,ICYCLE,CHOPT)
*
      ELSEIF (CHPATL.EQ.'LOCK') THEN
         CALL KUGETC(CHDIR,NCH)
         CALL RZLOCK(CHDIR)
*
      ELSEIF (CHPATL.EQ.'FREE') THEN
         CALL KUGETC(CHDIR,NCH)
         CALL RZFREE(CHDIR)
*
      ENDIF
  999 END
