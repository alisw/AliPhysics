*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXFCA(CHFUNC,NCH,JAD,IER)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *        To analize character string CHFUNC of length NCH        *
C.    *        CHFUNC may be the name of a COMIS function              *
C.    *        or a file name                                          *
C.    *                                                                *
C.    ******************************************************************
      COMMON/QUEST/IQUEST(100)
      CHARACTER*(*) CHFUNC
      CHARACTER*32 CHFILE
      INTEGER CSADDR
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
      JADF=0
      IQUEST(1)=0
#if defined(CERNLIB_COMIS)
      IF(NCH.LE.1.AND.CHFUNC(1:1).EQ.'0')GO TO 50
      IF(NCH.LE.0)GO TO 99
      IFILE=0
      DO 10 I=1,NCH-1
         IF(CHFUNC(I:I).EQ.'.')THEN
#endif
#if (defined(CERNLIB_COMIS))&&(defined(CERNLIB_IBM))
            CHFUNC(I:I)=' '
#endif
#if defined(CERNLIB_COMIS)
            IF(IFILE.EQ.0)IFILE=I-1
         ENDIF
  10  CONTINUE
  20  IF(IFILE.NE.0)THEN
         CHFILE='!FILE '//CHFUNC(1:NCH)
         CALL CSEXEC(CHFILE,IRET)
         DO 30 I=IFILE,1,-1
            IF(CHFUNC(I:I).EQ.'/'.OR.
     +         CHFUNC(I:I).EQ.BSLASH.OR.
     +         CHFUNC(I:I).EQ.']'.OR.
     +         CHFUNC(I:I).EQ.'~')THEN
                  I1=I+1
                  GO TO 40
            ENDIF
  30     CONTINUE
         I1=1
  40     CHFILE=CHFUNC(I1:IFILE)
         JADF=CSADDR(CHFILE)
         CHFUNC=CHFILE
         NCH=IFILE-I1+1
      ELSE
         JADF=CSADDR(CHFUNC)
      ENDIF
#endif
  50  JAD=JADF
*
  99  END
 
