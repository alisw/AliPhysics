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
      SUBROUTINE GXDZ
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      DZEBRA control commands                                   *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                  P.Zanarini  **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcxlun.inc"
      COMMON/GCLTMP/LTEMPA(3)
      COMMON/QUEST/IQUEST(100)
      CHARACTER*20 CHOPT
      CHARACTER*4 BNAME
      CHARACTER*32 CHPATL
      CHARACTER*255 DOCFIL
      CHARACTER*10  CHTOP
      CHARACTER*12  CHTDIR
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF(CHPATL.EQ.'SHOW')THEN
         CALL KUGETC(BNAME,NCH)
         CALL KUGETI(NUMBER)
         CALL KUGETC(CHOPT,NCH)
         LOC=LZLOC(IXSTOR,BNAME,NUMBER)
         IF(LOC.NE.0)THEN
            CALL DZSHOW('DZSHOW',IXSTOR,LOC,CHOPT,0,0,0,0)
         ENDIF
*
      ELSEIF(CHPATL.EQ.'SURV')THEN
         CALL KUGETC(BNAME,NCH)
         CALL KUGETI(NUMBER)
         LOC=LZLOC(IXSTOR,BNAME,NUMBER)
         IF(LOC.NE.0)THEN
            CALL DZSURV('DZSURV',IXSTOR,LOC)
         ENDIF
*
      ELSEIF(CHPATL.EQ.'SNAP')THEN
         CALL KUGETI(IDIV)
         CALL KUGETC(CHOPT,NCH)
         CALL DZSNAP('DZSNAP',IDIV,CHOPT)
*
      ELSEIF(CHPATL.EQ.'VERIFY')THEN
         CALL KUGETI(IDIV)
         CALL KUGETC(CHOPT,NCH)
         CALL DZVERI('DZVERI',IDIV,CHOPT)
*
      ELSEIF(CHPATL.EQ.'STORE')THEN
         CALL KUGETI(IXSTOR)
         CALL DZSTOR('DZSTOR',IXSTOR)
*
      ELSEIF (CHPATL.EQ.'DDIV') THEN
#if defined(CERNLIB_DZDOC)
         CALL KUGETI(IDIV)
         DOCFIL=' '
         CALL KUGETC(DOCFIL,NCH)
         NCH=MAX(NCH,1)
#endif
#if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DZDOC))
         CALL CUTOL(DOCFIL(1:NCH))
#endif
#if defined(CERNLIB_DZDOC)
         CHOPT='IN'
         IWDISP=1
         IWMETA=0
         ILOCNR=1
         CHTOP =' '
         CHTDIR=' '
         IF(DOCFIL.NE.' ') THEN
            CALL GXLUNF(10, 3, IFREE)
            IF(IFREE.NE.0) THEN
               LRECL = 256
               CALL RZOPEN(IFREE,CHTOP,DOCFIL(1:NCH),'W',LRECL,ISTAT)
               IF(ISTAT.NE.0) GOTO 999
               CALL RZFILE(IFREE,CHTOP,' ')
               LUNIT(IFREE) = 3
               CHTDIR='//'//CHTOP
            ENDIF
         ENDIF
*
*        Workstation ID is 1 for GEANT
*
         IWTYPE=IGIWTY(1)
         CALL DZDDIV(IDIV,LMAIN,CHTDIR,CHOPT,IWDISP,
     +               IWMETA,ILOCNR,IWTYPE)
         IF(DOCFIL.NE.' ') THEN
            CALL RZEND(CHTOP)
            CLOSE(IFREE)
            LUNIT(IFREE)=0
         ENDIF
#endif
#if !defined(CERNLIB_DZDOC)
*
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
10000    FORMAT(' *** GXDZ *** : Sorry, command not available',
     +          ' in this implementation')
#endif
      ELSEIF (CHPATL.EQ.'DISP') THEN
#if defined(CERNLIB_DZDOC)
         CALL KUGETC(BNAME,NCH)
         DOCFIL=' '
         CALL KUGETC(DOCFIL,NCH)
         NCH=MAX(NCH,1)
#endif
#if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DZDOC))
         CALL CUTOL(DOCFIL(1:NCH))
#endif
#if defined(CERNLIB_DZDOC)
         CALL KUGETI(NUMBER)
         CHOPT='N'
         IWDISP=1
         IWMETA=0
         ILOCNR=1
         CHTOP =' '
         CHTDIR=' '
         LOC=LZLOC(IXSTOR,BNAME,NUMBER)
         IF(LOC.NE.0)THEN
            IF(DOCFIL.NE.' ') THEN
               CALL GXLUNF(10, 3, IFREE)
               IF(IFREE.NE.0) THEN
                  LRECL = 256
                  CALL RZOPEN(IFREE,CHTOP,DOCFIL(1:NCH),
     +                        'W',LRECL,ISTAT)
                  IF(ISTAT.NE.0) GOTO 999
                  CALL RZFILE(IFREE,CHTOP,' ')
                  LUNIT(IFREE) = 3
                  CHTDIR='//'//CHTOP
               ENDIF
            ENDIF
*
*        Workstation ID is 1 for GEANT
*
            IWTYPE=IGIWTY(1)
            CALL DZDISP(IXSTOR,LOC,CHTDIR,CHOPT,IWDISP,
     +                  IWMETA,ILOCNR,IWTYPE)
            IF(DOCFIL.NE.' ') THEN
               CALL RZEND(CHTOP)
               CLOSE(IFREE)
               LUNIT(IFREE)=0
            ENDIF
         ENDIF
#endif
#if !defined(CERNLIB_DZDOC)
*
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
#endif
      ELSEIF (CHPATL.EQ.'DIRZ') THEN
#if defined(CERNLIB_DZDOC)
         DOCFIL=' '
         CALL KUGETC(DOCFIL,NCH)
         NCH=MAX(NCH,1)
#endif
#if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DZDOC))
         CALL CUTOL(DOCFIL(1:NCH))
#endif
#if defined(CERNLIB_DZDOC)
         CHOPT='N'
         IWDISP=1
         IWMETA=0
         ILOCNR=1
         IDZSTR=0
         CALL MZLINT(IDZSTR,'/GCLTMP/',LTEMPA,LTEMPA(3),LTEMPA)
         CALL DZDIRZ(IDZSTR,LTEMPA(3),1,DOCFIL(1:NCH),CHOPT,
     +               IWDISP,IWMETA,ILOCNR)
         LTEMPA(1)=0
#endif
#if !defined(CERNLIB_DZDOC)
         WRITE(CHMAIL,10000)
         CALL GMAIL(0,0)
#endif
      ENDIF
*
  999 END
