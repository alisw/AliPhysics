*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GPHOTI
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates cross-section in current material                  *
C     *  for PHOTOELECTRIC EFFECT using the empirical                  *
C     *  formula of L.Urban (see WRITE UP) for Z > 100                 *
C.    *                                                                *
C     *  otherwise Sandia parametrization                              *
C     *   F.Biggs & R. Lighthill, SAND87-0070, UC-34, May 1990         *
C     *                                                                *
C     *  NOTE:                                                         *
C     *  (1) Above 0.050 GeV ( 50 MeV) a cut is imposed as the         *
C     *   contribution from photoeffect is negligible for EGAM.GT.0.050*
C     *  (2) From 10KeV to 100KeV a special parametrization            *
C     *       has been introduced                                      *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author    L.Urban  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gctrak.inc"
#include "geant321/gconsp.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
C.
C.    ------------------------------------------------------------------
C.
C             universal energy cut-off
C
      SIG=0.
      EGAM=ELOW(IEKBIN)
      IF(EGAM.LE.0.05) THEN
C
         IF(Z.GT.100) THEN
C Below good old stuff
            IF(JMIXT.EQ.0)THEN
C
C             simple material (element)
C
               JPROB=LQ(JMA-4)
               SIG=GPHSIG(Z,EGAM)*AVO*DENS/A
            ELSE
C
C             compound/mixture
C
               NLMAT=Q(JMA+11)
               NLM=IABS(NLMAT)
               SIG=0.
               DO 10 I=1,NLM
                  J=JMIXT+I
                  AA=Q(J)
                  ZZ=Q(J+NLM)
                  WMAT=Q(J+2*NLM)
                  S=GPHSIG(ZZ,EGAM)
                  S=S*WMAT/AA
                  SIG=SIG+AVO*DENS*S
   10          CONTINUE
            ENDIF
C
C
         ELSEIF(Z.GE.1.) THEN
C New x-sections. Sandia parametrization
            JPROB=LQ(JMA-4)
            SIG = 0.0
            IF(EGAM.GE.1.E-8) SIG=GPHSG1(EGAM)*DENS
         ENDIF
      ENDIF
      IF(SIG.GT.0.)THEN
         Q(JPHOT+IEKBIN)=1./SIG
      ELSE
         Q(JPHOT+IEKBIN)=BIG
      ENDIF
C
      END
