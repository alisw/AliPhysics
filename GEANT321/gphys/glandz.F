*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/04 22/02/95  10.38.36  by  S.Giani
*-- Author :
      SUBROUTINE GLANDZ(Z,STEP,P,E,DEDX,DE,POTI,POTIL)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Energy straggling using a Monte Carlo model.                  *
C.    *  It can be used with or without delta ray generation.          *
C.    *                                                                *
C.    *   It is a NEW VERSION of the model , which reproduces          *
C.    *   the experimental data rather well.                           *
C.    *                                                                *
C.    *  Input : STEP  = current step-length (cm)                      *
C.    *                                                                *
C.    *  Output: DE    = actual energy loss (Gev)                      *
C.    *                 ( NOT the fluctuation DE/DX-<DE/DX> !)         *
C.    *                                                                *
C.    *     ==> Called by : GTELEC,GTHADR,GTMUON                       *
C.    *                                                                *
C.    *  Author      : L.Urban                                         *
C.    *  Date        : 28.04.1988       Last update :  1.02.90         *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gccuts.inc"
      PARAMETER (MAXRND=100)
      DIMENSION RNDM(MAXRND),APOIS(3),NPOIS(3),FPOIS(3)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION E1,E2,P2,B2,BG2,TM,DEC,W,WW,WWW,RR
      DOUBLE PRECISION X,AK,ALFA,EA,SA,AKL,REALK,FPOIS
      DOUBLE PRECISION ONE,HALF,ZERO
      PARAMETER (HMXINT=2**30)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (HMXINT=2**45)
#endif
      PARAMETER (ONE=1,HALF=ONE/2,ZERO=0)
      PARAMETER (RCD=0.40, RCD1=1.-RCD, PROBLM=0.01)
      PARAMETER( C1=4, C2=16)
*
*     ------------------------------------------------------------------
*
*     POTI=1.6E-8*Z**0.9
*     POTIL=LOG(POTI)
*
      IF(Z.GT.2.) THEN
         F2=2./Z
      ELSE
         F2=0.
      ENDIF
      F1=1.-F2
*
      E2 = 1.E-8*Z*Z
      E2L= LOG(E2)
      E1L= (POTIL-F2*E2L)/F1
      E1 = EXP(E1L)
*
*
      P2=P*P
      B2=P2/(E*E)
      BG2=P2/(AMASS*AMASS)
      IF(ITRTYP.EQ.2) THEN
         TM=P2/(E+EMASS)
         IF(CHARGE.LT.0.) TM=TM/2.
         TM=TM-POTI
         IF(TM.GT.DCUTE) TM=DCUTE
      ELSE
         TM=EMASS*P2/(0.5*AMASS*AMASS+EMASS*E)
         TM=TM-POTI
         IF(TM.GT.DCUTM) TM=DCUTM
      ENDIF
*
*
* *** Protection against negative TM    ---------------------
*     TM can be negative only for heavy particles with  a very
*     low kinetic energy (e.g. for proton with T  100-300 kev)
      TM=MAX(TM,ZERO)
*
      W  = TM+POTI
      WW = W/POTI
      WWW= 2.*EMASS*BG2
      WL = LOG(WWW)
      CSB=STEP*RCD1*DEDX/(WL-POTIL-B2)
      APOIS(1)=CSB*F1*(WL-E1L-B2)/E1
      APOIS(2)=CSB*F2*(WL-E2L-B2)/E2
*
      IF(TM.GT.0.) THEN
         APOIS(3)=RCD*DEDX*STEP*TM/(POTI*W*LOG(WW))
      ELSE
         APOIS(1)=APOIS(1)/RCD1
         APOIS(2)=APOIS(2)/RCD1
         APOIS(3)=0.
      ENDIF
*
*    calculate the probability of the zero energy loss
*
      APSUM=APOIS(1)+APOIS(2)+APOIS(3)
      IF(APSUM.LT.0.) THEN
         PROB=1.
      ELSEIF(APSUM.LT.50.) THEN
         PROB=EXP(-APSUM)
      ELSE
         PROB=0.
      ENDIF
*
*
*      do it differently if prob > problm  <====================
      IF(PROB.GT.PROBLM) THEN
         E0=1.E-8
         EMEAN=DEDX*STEP
         IF(TM.LE.0.) THEN
*      excitation only ....
            APOIS(1)=EMEAN/E0
*
            CALL GPOISS(APOIS,NPOIS,1)
            FPOIS(1)=NPOIS(1)
            DE=FPOIS(1)*E0
*
         ELSE
*         ionization only ....
            EM=TM+E0
            APOIS(1)=EMEAN*(EM-E0)/(EM*E0*LOG(EM/E0))
            CALL GPOISS(APOIS,NPOIS,1)
            NN=NPOIS(1)
            DE=0.
*
            IF(NN.GT.0) THEN
               RCORR=1.
               IF(NN.GT.MAXRND) THEN
                  RCORR=FLOAT(NN)/MAXRND
                  NN=MAXRND
*
               ENDIF
               W=(EM-E0)/EM
               CALL GRNDM(RNDM,NN)
               DO 10 I=1,NN
                  DE=DE+E0/(1.-W*RNDM(I))
   10          CONTINUE
               DE=RCORR*DE
*
            ENDIF
         ENDIF
         GOTO 999
      ENDIF
*
      IF(MAX(APOIS(1),APOIS(2),APOIS(3)).LT.HMXINT) THEN
         CALL GPOISS(APOIS,NPOIS,3)
         FPOIS(1)=NPOIS(1)
         FPOIS(2)=NPOIS(2)
         FPOIS(3)=NPOIS(3)
      ELSE
         DO 20 JPOIS=1, 3
            IF(APOIS(JPOIS).LT.HMXINT) THEN
               CALL GPOISS(APOIS(JPOIS),NPOIS(JPOIS),1)
               FPOIS(JPOIS)=NPOIS(JPOIS)
            ELSE
               CALL GRNDM(RNDM,2)
               FPOIS(JPOIS)=ABS(SQRT(-2.*LOG(RNDM(1)*ONE)
     +         *APOIS(JPOIS))*SIN(TWOPI*RNDM(2)*ONE)+APOIS(JPOIS))
            ENDIF
   20    CONTINUE
      ENDIF
 
*
*          Now we have all three numbers in REAL/DOUBLE
*          variables. REALK is actually an INTEGER that now may
*          exceed the machine representation limit for integers.
*
      DE=FPOIS(1)*E1+FPOIS(2)*E2
*
*     smear to avoid peaks in the energy loss (note: E1<<E2)
*
      IF(DE.GT.0.) THEN
         CALL GRNDM(RNDM,1)
         DE=DE+E1*(1.-2.*RNDM(1))
      ENDIF
*
      ALFA=1.
      REALK=0
      DEC=0.
*
      ANC=FPOIS(3)
      IF(ANC.GE.C2) THEN
         R=ANC/(C2+ANC)
         AN=ANC*R
         SN=C1*R
         CALL GRNDM(RNDM,2)
         RR=SQRT(-2.*LOG(RNDM(1)))
         PHI=TWOPI*RNDM(2)
         X=RR*COS(PHI)
         AK=AN+SN*X
         ALFA=WW*(C2+ANC)/(C2*WW+ANC)
         EA=AK*POTI*ALFA*LOG(ALFA)/(ALFA-1.)
         SA=SQRT(ABS(AK*ALFA*POTI*POTI-EA*EA/AK))
         AKL=(EA-C1*SA)/POTI
         IF(AK.LE.AKL) THEN
            X=RR*SIN(PHI)
            DEC=EA+SA*X
            REALK=AK+HALF-MOD(AK+HALF,ONE)
         ELSE
            ALFA=1.
         ENDIF
      ENDIF
      NN=NINT(FPOIS(3)-REALK)
      IF(NN.GT.MAXRND) THEN
         W=1.-ALFA/WW
         WW=POTI*ALFA
*
*     Here we take a gaussian distribution to avoid loosing
*     too much time in computing
*
         AVERAG=-LOG(1.-W)/W
         SIGMA =SQRT(NN*(1./(1.-W)-AVERAG**2))
         CALL GRNDM(RNDM,2)
         DEC   = DEC+WW*(NN*AVERAG+SIGMA*SQRT(-2.*LOG(RNDM(1)))*
     +           SIN(TWOPI*RNDM(2)))
         DE=DE+DEC
      ELSEIF(NN.GT.0) THEN
         W=1.-ALFA/WW
         WW=POTI*ALFA
         CALL GRNDM(RNDM,NN)
         DO 30 I=1,NN
            DEC=DEC+WW/(1.-W*RNDM(I))
   30    CONTINUE
         DE=DE+DEC
      ENDIF
*
  999 END
