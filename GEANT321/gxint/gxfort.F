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
      SUBROUTINE GXFORT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      COMIS/FORTRAN control commands                            *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcomis.inc"
#include "geant321/gcxlun.inc"
#include "geant321/gcunit.inc"
*
      INTEGER CSADDR
      CHARACTER*80 CHFILE
      CHARACTER*32 CHPATL
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF(CHPATL.EQ.'FORTRAN')THEN
#if defined(CERNLIB_COMIS)
         CALL KUGETC(CHFILE,NCH)
         IF(NCH.EQ.1.AND.CHFILE(1:1).EQ.'-')THEN
            JUINIT=0
            JUGEOM=0
            JUKINE=0
            JUSTEP=0
            JUOUT =0
            JULAST=0
         ELSEIF(NCH.GT.0) THEN
            CALL CLTOU(CHFILE)
            CALL GXFCA(CHFILE,NCH,JAD,1)
            JUINIT=CSADDR('UGINIT')
            JUGEOM=CSADDR('UGEOM')
            JUKINE=CSADDR('GUKINE')
            JUSTEP=CSADDR('GUSTEP')
            JUOUT =CSADDR('GUOUT')
            JULAST=CSADDR('UGLAST')
         ENDIF
*
#endif
      ENDIF
*
      END
