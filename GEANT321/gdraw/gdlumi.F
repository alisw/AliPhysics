*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:30  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GDLUMI(LINCOL,APROSC)
*
******************************************************************
*                                                                *
* RAY-TRACING                                                    *
* Function: Recieve light intensity and angle with surface       *
*           Compute colour shading for the surface               *
*                                                                *
* I/O parameters:                                                *
* LINCOL = colour code                                           *
* APROSC = inclination between light rays and surface            *
*                                                                *
* Author: S. Giani 1992                                          *
*                                                                *
******************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
*
#include "geant321/gctrak.inc"
#include "geant321/gcpixe.inc"
#include "geant321/gcrayt.inc"
 
      DATA IAAA /456789/
      SAVE IAAA
*
      IF(LINCOL.EQ.2)THEN
         LINCOL=16
      ELSEIF(LINCOL.EQ.3)THEN
         LINCOL=66
      ELSEIF(LINCOL.EQ.4)THEN
         LINCOL=116
      ELSEIF(LINCOL.EQ.5)THEN
         LINCOL=41
      ELSEIF(LINCOL.EQ.6)THEN
         LINCOL=141
      ELSEIF(LINCOL.EQ.7)THEN
         LINCOL=91
      ELSE
         print *,'Warning, color not supported'
         GOTO 998
      ENDIF
*
      LINSAV=LINCOL
      IAAA=MOD(IAAA*1237,1000000)+1
      if(iaaa.lt.1.or.iaaa.gt.1000000)print *,'iaaa<0'
      QIAAA=IAAA
*
      SMIN=-0.001
      SMAX=0.04
      DO 1 I=1,25
       IF(APROSC.GT.SMIN.AND.APROSC.LE.SMAX)THEN
        QD=QIAAA/1000000.
        QCD=QD*(SMAX-SMIN)
        IF((APROSC-SMIN).GT.QCD)THEN
           LINCOL=LINCOL+I
        ELSE
           IF(I.GT.1)THEN
            LINCOL=LINCOL+I-1
           ELSE
            LINCOL=LINCOL+I
           ENDIF
        ENDIF
         GOTO 998
       ELSE
         SMIN=SMAX
         SMAX=SMAX+0.04
       ENDIF
 1    CONTINUE
*
      print *,aprosc,'=aprosc'
*
 998  CONTINUE
*      CALL ISFACI(LINCOL)
*      CALL ISPLCI(LINCOL)
**      IF(ISWIT(9).NE.54321)CALL ISPMCI(LINCOL)
       CALL ISPMCI(LINCOL)
       ICOLOR=LINCOL
*
 999  END
 
