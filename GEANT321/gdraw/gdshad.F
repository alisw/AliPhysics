*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
#if defined(CERNLIB_CG)
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDSHAD(LINCOL,APROSC)
*
******************************************************************
*                                                                *
* Function: Recieve light intensity and angle with surface       *
*           Compute colour shading for the surface               *
*                                                                *
* I/O parameters:                                                *
* LINCOL = colour code                                           *
* APROSC = inclination between light rays and surface            *
*                                                                *
* Author: S. Giani                                               *
*                                                                *
******************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/cggpar.inc"
#include "geant321/cghpar.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
*
#include "geant321/gctrak.inc"
#include "geant321/gcrayt.inc"
*     DATA AAA /1./
*     SAVE AAA
*
*      print *,aprosc
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
*        print *,'color not supported'
         GOTO 20
      ENDIF
*
**      AAA=AAA+1
*******      AAA=AAA*123.456789
*******      AAA=AINT(AAA)
*******      QD=AAA
*
      SMIN=-1.0001
      SMAX=-1.+0.08
      DO 10 I=1,25
         IF(SMAX.GT.0.99)SMAX=1.0001
         IF(APROSC.GT.SMIN.AND.APROSC.LE.SMAX)THEN
**        QD=RNDM(AAA)
**        QCD=QD*(SMAX-SMIN)
**        IF((APRO-SMIN).GT.QCD)THEN
            LINCOL=LINCOL+(I+1)/1.5
**        ELSE
**           IF(I.GT.1)THEN
**            LINCOL=LINCOL+I-1
**           ELSE
**            LINCOL=LINCOL+I
**           ENDIF
**        ENDIF
            GOTO 20
         ELSE
            SMIN=SMAX
            SMAX=SMAX+0.08
         ENDIF
   10 CONTINUE
*
*     print *,aprosc,'=aprosc'
*
   20 CONTINUE
      CALL ISFACI(LINCOL)
      CALL ISPLCI(LINCOL)
      CALL ISPMCI(LINCOL)
*
  999 END
#endif
