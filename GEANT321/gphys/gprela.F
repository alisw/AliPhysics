*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:32  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.23  by  S.Giani
*-- Author :
      SUBROUTINE GPRELA
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Initialise energy loss due to direct pair-production  and     *
C.    *  nuclear interactions by muons.                                *
C.    *                                                                *
C     *       A, Direct pair production  <=======                      *
C     *                                                                *
C.    *  The DE/DX expression of MANDO and RONCHI(NUOVO CIMENTO        *
C.    *  9(1952),517) is used which attempts to account for the effect *
C.    *  of atomic electron screening. The screening correction is     *
C.    *  applied above the energy limit defined by C. RICHARD-SERRE    *
C.    *  (CERN 71-18).                                                 *
C.    *                                                                *
C     *       B, Nuclear interactions  <=======                        *
C.    *                                                                *
C.    *  the following expression derived from the Williams-Weizsacker *
C.    *  relation for the virtual photon flux is used :                *
C.    *                                                                *
C.    *      DE/DX = FACTOR * XSEC * E                                 *
C.    *                                                                *
C.    *      where,                                                    *
C.    *      FACTOR  =  2. * N * ALPHA / PI                            *
C.    *      XSEC    =  total photo-nuclear cross-section assumed to   *
C.    *                 be independent of energy. Value of 140 mubarns *
C.    *                 calculated by averaging the results of Hesse   *
C.    *                 et al(Phys. Rev. Lett 25(1970),613).           *
C.    *                                                                *
C.    *    ==>Called by : GPHYSI                                       *
C.    *       Author    G.Patrick  *********                           *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gcmulo.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"
#include "geant321/gctrak.inc"
      DATA XSEC  /140.E-30/
      DATA FACTOR/2.7976238E+21/
C.
C.    ------------------------------------------------------------------
C.
      IF(Z.LT.1.) GOTO 999
      ICHAN=IEKBIN
      T    = ELOW(ICHAN)
      T1   = 10.**(EKBIN(1)+(ICHAN-0.5)/GEKA)
      IF(T.LT.1.)GO TO 999
      E    = T+EMMU
      E1   = T1+EMMU
      IF(JMIXT.EQ.0)THEN
C
C           Element
C
         DEDX = GPRELM(Z,T,PPCUTM)
         DEDX = AVO*DENS*DEDX/A
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            DEDX1 = GPRELM(Z,T1,PPCUTM)
            DEDX1 = AVO*DENS*DEDX1/A
         ENDIF
      ELSE
C
C          Compound/Mixture
C
         NLMAT = Q(JMA+11)
         NLM   = IABS(NLMAT)
         DEDX  = 0.
         DEDX1 = 0.
         DO 10 L=1,NLM
            J    = JMIXT+NLM+L
            AA   = Q(J-NLM)
            ZZ   = Q(J)
            WMAT = Q(J+NLM)
            S    = GPRELM(ZZ,T,PPCUTM)
            S    = WMAT*S/AA
            DEDX = DEDX+AVO*DENS*S
*
* *** auxiliary integration point for Range tables
            IF(ICHAN.NE.NEK1) THEN
               S = GPRELM(ZZ,T1,PPCUTM)
               S = WMAT*S/AA
               DEDX1 = DEDX1+AVO*DENS*S
            ENDIF
   10    CONTINUE
      ENDIF
C
C             Nuclear interactions
C
      IF(IMUNU.EQ.0.AND.E.GE.10.) THEN
         DENU=DENS*FACTOR*XSEC*E
*
* *** auxiliary integration point for Range tables
         IF(ICHAN.NE.NEK1) THEN
            DENU1=DENS*FACTOR*XSEC*E1
         ENDIF
      ELSE
         DENU=0.0
         DENU1=0.0
      ENDIF
      IF(DEDX.LT.0.)DEDX=0.
      IF(DENU.LT.0.)DENU=0.
      JEL2=LQ(JMA-2)
      Q(JEL2+ICHAN)=Q(JEL2+ICHAN)+DEDX+DENU
*
* *** auxiliary integration point for Range tables
      IF(ICHAN.NE.NEK1) THEN
         IF(DEDX1.LT.0.)DEDX1=0.
         IF(DENU1.LT.0.)DENU1=0.
         WS(NEKBIN*2+ICHAN)=WS(NEKBIN*2+ICHAN)+DEDX1+DENU1
      ENDIF
C
  999 END
