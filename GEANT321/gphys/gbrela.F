*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GBRELA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Initialize the bremsstrahlung energy loss tables              *
C.    *            ( See L.Urban report)                               *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author    G.Patrick, L.Urban  *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gctrak.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcmate.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcphys.inc"
      DIMENSION DEM(8)
      DATA DEM/-6.76228E-1,4.35611E-1,-4.69224E-1,-7.05133E-3,
     +        1.31200,-8.07563E-1,-6.96166E-1,3.62549E-1/
C.
C.    ------------------------------------------------------------------
C.
C            Number of constituents(ie. element,mixture or compound).
C
      NLMAT = Q(JMA+11)
      NLM    = IABS(NLMAT)
      IF (NLMAT.EQ.0)                             GO TO 999
      IF(Z.LT.1.)                                 GO TO 999
C
      JEL1=LQ(JMA-1)
      JEL2=LQ(JMA-2)
      ICHAN=IEKBIN
      XE     = ELOW(ICHAN)
      XE1    = 10.**(EKBIN(1)+(ICHAN-0.5)/GEKA)
C
C==========>  A,Bremsstrahlung losses for electrons
C             special treatment for compounds/mixtures
C
      IF(JMIXT.EQ.0)THEN
C
C             simple material (element)
C
         DEDX=GBRELE(Z,XE,BCUTE)
         DEDX=AVO*DENS*DEDX/A
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            DEDX1=GBRELE(Z,XE1,BCUTE)
            DEDX1=AVO*DENS*DEDX1/A
         ENDIF
      ELSE
C
C             compound/mixture
C
         DEDX=0.
         DEDX1=0.
         DO 10 L=1,NLM
            J = JMIXT+NLM+L
            AA = Q(J-NLM)
            ZZ = Q(J)
            WMAT= Q(J+NLM)
            S = GBRELE(ZZ,XE,BCUTE)
            S = WMAT*S/AA
            DEDX= DEDX+AVO*DENS*S
*
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               S = GBRELE(ZZ,XE1,BCUTE)
               S = WMAT*S/AA
               DEDX1= DEDX1+AVO*DENS*S
            ENDIF
   10    CONTINUE
      ENDIF
C
      IF(DEDX.LT.0.)DEDX=0.
      Q(JEL1+ICHAN)=Q(JEL1+ICHAN)+DEDX
      FACT = GBFLOS(XE,BCUTE)
      Q(JEL1+ICHAN+NEK1)=Q(JEL1+ICHAN+NEK1)+FACT*DEDX
*
      IF(ISTRA.NE.0) Q(JEL1+ICHAN+2*NEK1) = DEDX
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         IF(DEDX1.LT.0.)DEDX1=0.
         WS(ICHAN)=WS(ICHAN)+DEDX1
         FACT = GBFLOS(XE1,BCUTE)
         WS(NEKBIN+ICHAN)=WS(NEKBIN+ICHAN)+FACT*DEDX1
      ENDIF
C
C===========>  B, Bremsstrahlung losses for muons
C
      IF(XE.GE.1.)GO TO 20
      CONS1=3.26485E-8*Z*Z*DENS/A
      CONS2=LOG(2.70511)-LOG(A)/3.
      R=EMMU/(XE+EMMU)
      RL=LOG(R)
      DEDX=DEM(1)+RL*DEM(3)+R*(DEM(5)+R*DEM(7))
      DEDX=DEDX+CONS2*(DEM(2)+RL*DEM(4)+R*(DEM(6)+R*DEM(8)))
      DEDX=CONS1*XE*DEDX
      IF(DEDX.LT.0.)DEDX=0.
      Q(JEL2+ICHAN)=Q(JEL2+ICHAN)+DEDX
*
      IF(ISTRA.NE.0) Q(JEL2+ICHAN+NEK1) = DEDX
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         R=EMMU/(XE1+EMMU)
         RL=LOG(R)
         DEDX1=DEM(1)+RL*DEM(3)+R*(DEM(5)+R*DEM(7))
         DEDX1=DEDX1+CONS2*(DEM(2)+RL*DEM(4)+R*(DEM(6)+R*DEM(8)))
         DEDX1=CONS1*XE1*DEDX1
         IF(DEDX1.LT.0.)DEDX1=0.
         WS(NEKBIN*2+ICHAN) = WS(NEKBIN*2+ICHAN)+DEDX1
      ENDIF
      GO TO 999
C
   20 IF(JMIXT.EQ.0)THEN
C
C             simple material (element)
C
         DEDX=GBRELM(Z,XE,BCUTM)
         DEDX=AVO*DENS*DEDX/A
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            DEDX1=GBRELM(Z,XE1,BCUTM)
            DEDX1=AVO*DENS*DEDX1/A
         ENDIF
      ELSE
C
C             compound/mixture
C
         DEDX=0.
         DEDX1=0.
         DO 30 L=1,NLM
            J = JMIXT+NLM+L
            AA = Q(J-NLM)
            ZZ = Q(J)
            WMAT= Q(J+NLM)
            S = GBRELM(ZZ,XE,BCUTM)
            S = WMAT*S/AA
            DEDX= DEDX+AVO*DENS*S
*
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               S = GBRELM(ZZ,XE1,BCUTM)
               S = WMAT*S/AA
               DEDX1= DEDX1+AVO*DENS*S
            ENDIF
   30    CONTINUE
      ENDIF
C
      IF(DEDX.LT.0.)DEDX=0.
      Q(JEL2+ICHAN)=Q(JEL2+ICHAN)+DEDX
*
      IF(ISTRA.NE.0) Q(JEL2+ICHAN+NEK1) = DEDX
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         IF(DEDX1.LT.0.)DEDX1=0.
         WS(NEKBIN*2+ICHAN)=WS(NEKBIN*2+ICHAN)+DEDX1
      ENDIF
C
  999 END
