*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:33  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GRANGI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the stopping range , in cm.                        *
C.    *     (i.e. the maximum step due to the continuous energy loss)  *
C.    *      The stopping range is the integral of the inverse of the  *
C.    *      DE/DX table.                                              *
C.    *     There are 4 tables : electron,positron,muon,proton         *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author      M.Maire   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION RANGMU, RANGEL, RANGPO, RANGPR
      DOUBLE PRECISION DDXEL1, DDXMU1, DDXPR1, DDXPO1
      DOUBLE PRECISION DDXEL2, DDXMU2, DDXPR2, DDXPO2
      DOUBLE PRECISION DDXEL3, DDXMU3, DDXPR3, DDXPO3
      DOUBLE PRECISION HFACT, ONE, TEN, HLOG10, HFACT2, HFACT3
#endif
      PARAMETER (ONE=1,TEN=10)
*
      HLOG10 = LOG(TEN)
*
      JRANEL = LQ(JMA-15)
      JRANPO = LQ(JMA-15)+NEK1
      JRANMU = LQ(JMA-16)
      JRANPR = LQ(JMA-16)+NEK1
*
      IF(Z.LT.1) THEN
         DO 10  IEKBIN=1,NEK1
            Q(JRANEL+IEKBIN)=BIG
            Q(JRANPO+IEKBIN)=BIG
            Q(JRANMU+IEKBIN)=BIG
            Q(JRANPR+IEKBIN)=BIG
   10    CONTINUE
*
      ELSE
         SMALL = 1./BIG
         JLOSEL = LQ(JMA-1)
         JLOSPO = LQ(JMA-1)+NEK1
         JLOSMU = LQ(JMA-2)
         JLOSPR = LQ(JMA-3)
*
         JWSPEL =         -1
         JWSPPO = NEKBIN  -1
         JWSPMU = NEKBIN*2-1
         JWSPPR = NEKBIN*3-1
*
         RANGPO = 0.
         RANGMU = 0.
         RANGEL = 0.
         RANGPR = 0.
*
         Q(JRANEL+1)=0.
         Q(JRANPO+1)=0.
         Q(JRANMU+1)=0.
         Q(JRANPR+1)=0.
*
         HFACT = ONE/(6.*GEKA)
*
         DDXEL3 = ELOW(1)*HLOG10/MAX(Q(JLOSEL+1),SMALL)
         DDXMU3 = ELOW(1)*HLOG10/MAX(Q(JLOSMU+1),SMALL)
         DDXPR3 = ELOW(1)*HLOG10/MAX(Q(JLOSPR+1),SMALL)
         DDXPO3 = ELOW(1)*HLOG10/MAX(Q(JLOSPO+1),SMALL)
*
         HFACT2 = HLOG10*TEN**(EKBIN(1)+0.5/GEKA)
*
         DO 20 IEKBIN=2,NEK1
*
            HFACT3 = ELOW(IEKBIN)*HLOG10
*
* *** Electrons
*
            DDXEL1 = DDXEL3
            DDXEL3 = HFACT3/MAX(Q(JLOSEL+IEKBIN),SMALL)
            DDXEL2 = HFACT2/MAX(WS(JWSPEL+IEKBIN),SMALL)
            RANGEL = RANGEL+DDXEL1+4.*DDXEL2+DDXEL3
            Q(JRANEL+IEKBIN)=RANGEL*HFACT
*
* *** Positons
*
            DDXPO1 = DDXPO3
            DDXPO3 = HFACT3/MAX(Q(JLOSPO+IEKBIN),SMALL)
            DDXPO2 = HFACT2/MAX(WS(JWSPPO+IEKBIN),SMALL)
            RANGPO = RANGPO+DDXPO1+4.*DDXPO2+DDXPO3
            Q(JRANPO+IEKBIN)=RANGPO*HFACT
*
* *** Muons
*
            DDXMU1 = DDXMU3
            DDXMU3 = HFACT3/MAX(Q(JLOSMU+IEKBIN),SMALL)
            DDXMU2 = HFACT2/MAX(WS(JWSPMU+IEKBIN),SMALL)
            RANGMU = RANGMU+DDXMU1+4.*DDXMU2+DDXMU3
            Q(JRANMU+IEKBIN)=RANGMU*HFACT
*
* *** Protons
*
            DDXPR1 = DDXPR3
            DDXPR3 = HFACT3/MAX(Q(JLOSPR+IEKBIN),SMALL)
            DDXPR2 = HFACT2/MAX(WS(JWSPPR+IEKBIN),SMALL)
            RANGPR = RANGPR+DDXPR1+4.*DDXPR2+DDXPR3
            Q(JRANPR+IEKBIN)=RANGPR*HFACT
*
            HFACT2 = HLOG10*TEN**(EKBIN(1)+(IEKBIN-0.5)/GEKA)
*
   20    CONTINUE
      ENDIF
*
      END
