*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.38  by  S.Giani
*-- Author :
      FUNCTION GMUSIG(E,E1,COST)
C
C *** CALCULATION OF MUON-NUCLEUS CROSS SECTIONS IN CURRENT MATERIAL ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C CALLED BY : GMUNUI
C ORIGIN : F.CARMINATI, H.FESEFELDT
C
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
C
      REAL NU,K,MMU
      PARAMETER (ALPHA=0.007299)
C
C
C --- GET MUON MASS FROM GEANT PARTICLE DATA BANK ---
      LNVE=LQ(JPART-5)
      MMU=Q(LNVE+7)
C
      GMUSIG=0.
      IF(E1.LE.MMU) GO TO 9999
      IF(COST.GE.1.) THEN
         GMUSIG=BIG
         GOTO 9999
      ELSE
         TAGTB2 = (1.-COST)/SQRT((1.-COST)*(1.+COST))
      ENDIF
      P      = E*E-MMU*MMU
      P1     = E1*E1-MMU*MMU
      IF(P.LT.0..OR.P1.LT.0.) GOTO 9999
      P      = SQRT(P)
      P1     = SQRT(P1)
      Q2     = 2.*(E*E1-P*P1*COST-MMU*MMU)
      IF(Q2.LE.0.) GOTO 9999
      W2 = PMASS*( PMASS + 2.0*(E-E1) ) - Q2
      IF ( W2.LE.0.0 ) GOTO 9999
      NU     = E-E1
      K      = NU+0.5*Q2/NU
      EPS    = 1./(1.+2.*TAGTB2*TAGTB2*(Q2+NU*NU)/Q2)
      IF(EPS.GE.1) THEN
         GMUSIG=BIG
      ELSE
         GAM    = ALPHA*K*E1/(PI*Q2*E*(1.-EPS))
         ST     = 0.12
         SL     = 0.3*(1.-Q2/(NU*1.868))*ST
         IF(SL.LT.0.) SL=0.
         SLST   = ST+EPS*SL
         GMUSIG = MAX(GAM*SLST,0.0)
      ENDIF
C
 9999 CONTINUE
      END
