*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRELE(EEL,CHARGE,JMA,DEDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Compute ion losses for electron/positron                 *
C.    *                                                                *
C.    *    ==>Called by : GDRELA                                       *
C.    *       Author    G.Patrick *********                            *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"
C
      DATA CONS/0.153536E-3/
C.
C.    ------------------------------------------------------------------
C.
      JPROB=LQ(JMA-4)
C
      GAM=EEL/EMASS + 1.
      GAM2=GAM*GAM
      T=GAM-1.
      DEDX=0.
      IF(T.LE.0.)GO TO 99
      TCME=DCUTE/EMASS
      BET2=1.-1./GAM2
C     ------------------------------
      IF(CHARGE.GT.0.) THEN
         Y=1./(1.+GAM)
         D=TCME
         IF(T.LT.TCME) D=T
         D2=D*D/2.
         D3=2.*D2*D/3.
         D4=D2*D2
         F=LOG(T*D)-BET2*(T+2.*D-Y*(3.*D2
     *    +Y*(D-D3+Y*(D2-T*D3+D4))))/T
C
      ELSE
        D=TCME
         IF(T.LT.2.*TCME) D=0.5*T
         F=-1.-BET2+LOG((T-D)*D)+T/(T-D)
     *    +(0.5*D*D+(1.+2.*T)*LOG(1.-D/T))/GAM2
      ENDIF
C
      POTL=Q(JPROB+11)
      FAC=Q(JPROB+12)
      C=Q(JPROB+13)
      X0=Q(JPROB+14)
      X1=Q(JPROB+15)
      AA=Q(JPROB+16)
C
      X=LOG(GAM2-1.)/4.606
      DEL=0.
      IF(X.GT.X0)THEN
         DEL=4.606*X+C
         IF(X.LE.X1)DEL=DEL+AA*(X1-X)**3.
      ENDIF
C
      DEDX=CONS*FAC*(LOG(2.*T+4.)-2.*POTL+F-DEL)/BET2
      IF(DEDX.LT.0.)DEDX=0.
C
  99  RETURN
      END
