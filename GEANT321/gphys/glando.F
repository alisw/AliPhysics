*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GLANDO(IMODE,STEP,Z,A,RHO,P,E,XMASS,DE,IFLAG)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *                                                                *
C.    *  GLANDO                                                        *
C.    *  ------                                                        *
C.    *                                                                *
C.    *  Energy straggling using Gaussian, Landau & Vavilov theories.  *
C.    *                                                                *
C.    *  Input                                                         *
C.    *  -----                                                         *
C.    *  IMODE  =  3   Landau(RANLAN)  sampling                        *
C.    *         =  4   Landau(GENLAN)  sampling                        *
C.    *         =  5   Vavilov(DINVAV) sampling                        *
C.    *         =  6   Gaussian sampling                               *
C.    *         =  2   Automatic selection of relevant distribution    *
C.    *                                                                *
C.    *  STEP   =  current step-length (cm)                            *
C.    *                                                                *
C.    *  Output                                                        *
C.    *  ------                                                        *
C.    *  DE     =  DE/DX - <DE/DX>     (GeV)                           *
C.    *                                                                *
C.    *  IFLAG  =  3   Landau(RANLAN) sampling used                    *
C.    *         =  4   Landau(GENLAN) sampling used                    *
C.    *         =  5   Vavilov(GVAVIV) sampling used                   *
C.    *         =  6   Gaussian sampling used                          *
C.    *                                                                *
C.    *  Warning                                                       *
C.    *  -------                                                       *
C.    *  Only Landau sampling should be used since this has been well  *
C.    *  tested whereas both Vavilov and Gaussian sampling are being   *
C.    *  developed.                                                    *
C.    *                                                                *
C.    *  Author      : G.N. Patrick                                    *
C.    *  Date        : 03.05.1985                                      *
C.    *  Last update : 09.09.1985                                      *
C.    *                                                                *
C.    ******************************************************************
C.
      PARAMETER (EULER=0.577215)
      PARAMETER (P1=.60715,P2=.67794,P3=.52382E-1,P4=.94753,
     +           P5=.74442,P6=1.1934)
#include "geant321/gconsp.inc"
      DIMENSION RNDM(2)
      FLAND(X) = P1+P6*X+(P2+P3*X)*EXP(P4*X+P5)
C.
C.    ------------------------------------------------------------------
C.
      IF(STEP.LT.1.E-7)THEN
         DE=0.
         IFLAG=0
         RETURN
      ENDIF
C
C     Calculate xi factor (KeV).
C
      BETA   = P/E
      GAMMA  = E/XMASS
      XI     = (153.5*Z*STEP*RHO)/(A*BETA*BETA)
C
C     Maximum energy transfer to atomic electron (KeV).
C
      ETA    = BETA*GAMMA
      ETASQ  = ETA*ETA
      RATIO  = EMASS/XMASS
      F1     = 2.*EMASS*ETASQ
      F2     = 1.+2.*RATIO*GAMMA+RATIO*RATIO
      EMAX   = F1*1.E+6/F2
C
C     Calculate Kappa significance ratio.
C
      CAPPA  = XI/EMAX
C
C     Choose correct function if IMODE set to automatic selection.
C
      IMODEI = IMODE
      IF (IMODEI.EQ.2)    THEN
         IF (CAPPA.LT.0.01)                              IMODEI = 4
*****    IF (CAPPA.GE.0.01.AND.CAPPA.LE.10.)             IMODEI = 3
         IF (CAPPA.GE.0.01.AND.CAPPA.LE.10.)             IMODEI = 5
         IF (CAPPA.GT.10.)                               IMODEI = 6
      ENDIF
C
C     +---------------------------------------------------------------+
C     I Sample lambda variable from Kolbig/Schorr Landau distribution I
C     +---------------------------------------------------------------+
C
   10 CONTINUE
      IF (IMODEI.EQ.3) THEN
         XMEAN = -BETA**2-LOG(XI/EMAX)+EULER-1.
         XLAMX = FLAND(XMEAN)
   20    CALL GRNDM(RNDM,1)
         IF( RNDM(1) .GT. 0.980 ) GO TO 20
         XLAMB = GLANDR(RNDM(1))
         IF(XLAMB.GT.XLAMX) GO TO 20
         IFLAG  = 3
C
C            +----------------------------------------------+
C            I  Sample lambda variable from James & Hancock I
C            I  Landau distribution                         I
C            +----------------------------------------------+
C
 
      ELSEIF (IMODEI.EQ.4) THEN
         XMEAN = -BETA**2-LOG(XI/EMAX)+EULER-1.
         XLAMX = FLAND(XMEAN)
  25     CALL GLANDG(XLAMB)
         IF(XLAMB.GT.XLAMX) GO TO 25
         IFLAG = 4
C
C            +---------------------------------------------------------+
C            I Sample lambda variable (Landau not Vavilov) from        I
C            I Rotondi&Montagna&Kolbig Vavilov distribution            I
C            +---------------------------------------------------------+
C
      ELSEIF (IMODEI.EQ.5) THEN
C
C            Keep within computation bounds.
C
         RKA = CAPPA
         BE2 = BETA*BETA
         CALL GRNDM(RNDM,1)
         XLAMB = GVAVIV(RKA,BE2,RNDM(1))
         IFLAG = 5
C
C     +-----------------------------------+
C     I Sample from Gaussian distribution I
C     +-----------------------------------+
C
      ELSEIF (IMODEI.EQ.6) THEN
         SIGMA  = XI*EMAX*(1.-(BETA*BETA/2.))
         SIGMA  = SQRT(SIGMA)
   30    CALL GRNDM(RNDM,2)
         IF(RNDM(1).LE.0.)GO TO 30
         F1     = -2.*LOG(RNDM(1))
         DEKEV  = SIGMA*SQRT(F1)*COS(2.*PI*RNDM(2))
         IFLAG  = 6
         GOTO 40
      ENDIF
C
C     Calculate DE/DX - <DE/DX>
C
      DEKEV  = XI*(XLAMB+BETA*BETA+LOG(XI/EMAX)-EULER+1.)
C
   40 DE=DEKEV*1.E-6
      END
