*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
*
      SUBROUTINE CGFAC2(TVVX,TVVY,TVVZ,IN)
*
********************************************************************
*                                                                  *
* Function: Order edges in world coordinates for every face        *
*                                                                  *
* I/O parameters:                                                  *
* TVVX,TVVY,TVVZ = set of edge's coordinates                       *
* IN = number of edges *2.                                         *
*                                                                  *
* Author: S. Giani                                                 *
*                                                                  *
********************************************************************
*
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cghpar.inc"
#include "geant321/cgctra.inc"
#include "geant321/cgcedg.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcspee.inc"
*
      DIMENSION VVX(500),VVY(500),VVZ(500)
      DIMENSION TVVX(500),TVVY(500),TVVZ(500)
      DIMENSION IZ(500)
*
      JSC=0
      VVX(1)=TVVX(1)
      VVY(1)=TVVY(1)
      VVZ(1)=TVVZ(1)
      VVX(2)=TVVX(2)
      VVY(2)=TVVY(2)
      VVZ(2)=TVVZ(2)
      IZ(1)=1
      IZ(2)=2
      KZ=3
      DO 143 II=3,IN,2
       DO 144 JJ=3,IN
        DO 9 KK=1,KZ-1
         IF(JJ.EQ.IZ(KK))GOTO 144
 9      CONTINUE
          C11=ABS(TVVX(JJ)-VVX(II-1))
          C12=ABS(TVVY(JJ)-VVY(II-1))
          C13=ABS(TVVZ(JJ)-VVZ(II-1))
          IF(C11.LT..001.AND.C12.LT..001.AND.C13.LT..001)THEN
            VVX(II)=TVVX(JJ)
            VVY(II)=TVVY(JJ)
            VVZ(II)=TVVZ(JJ)
            IZ(KZ)=JJ
            AJ=JJ*.5
            IAJ=AJ
            PDAJ=AJ-IAJ
            IF(PDAJ.GT..01)THEN
              VVX(II+1)=TVVX(JJ+1)
              VVY(II+1)=TVVY(JJ+1)
              VVZ(II+1)=TVVZ(JJ+1)
              IZ(KZ+1)=JJ+1
              KZ=KZ+2
            ELSE
              VVX(II+1)=TVVX(JJ-1)
              VVY(II+1)=TVVY(JJ-1)
              VVZ(II+1)=TVVZ(JJ-1)
              IZ(KZ+1)=JJ-1
              KZ=KZ+2
            ENDIF
            GOTO 143
          ENDIF
 144   CONTINUE
       DO 11 JJA=3,IN
        DO 10 I=1,KZ-1
        IF(JJA.EQ.IZ(I))GOTO 11
 10     CONTINUE
        JJAO=JJA
        GOTO 12
 11     CONTINUE
        PRINT *,'Error in CGFACO'
 12     CONTINUE
        VVX(II)=TVVX(JJAO)
        VVY(II)=TVVY(JJAO)
        VVZ(II)=TVVZ(JJAO)
        VVX(II+1)=TVVX(JJAO+1)
        VVY(II+1)=TVVY(JJAO+1)
        VVZ(II+1)=TVVZ(JJAO+1)
        IZ(KZ)=JJAO
        IZ(KZ+1)=JJAO+1
        KZ=KZ+2
        JSC=II-1
 143  CONTINUE
      DO 145 KK=1,IN
       TVVX(KK)=VVX(KK)
       TVVY(KK)=VVY(KK)
       TVVZ(KK)=VVZ(KK)
 145  CONTINUE
*
 999  END
 
