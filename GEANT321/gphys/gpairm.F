*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.22  by  S.Giani
*-- Author :
      SUBROUTINE GPAIRM
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates direct pair production by muons                     *
C.    *                                                                *
C.    *    ==>Called by : GTMUON                                       *
C.    *       Author    L.Urban  *********                             *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcmate.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gcphys.inc"
#include "geant321/gccuts.inc"
      DIMENSION RNDM(2)
      LOGICAL ROTATE
C
C         CVM= 3*SQRT(e)*EMMU/4
C         EM6= 6*EMMU**2
      DATA CVM,EM6/0.130652,0.066983/
      DATA AL10T/9.212/
C.
C.    ------------------------------------------------------------------
C.
      IF(GEKIN.LE.PPCUTM)GO TO 900
      EEM1=GETOT
      KCASE=NAMEC(6)
C
      VMIN=4*EMASS/EEM1
      VMAX=1.-CVM*Z**0.333333/EEM1
      IF(VMAX.LE.VMIN)GO TO 900
      VC  = PPCUTM/EEM1
      ALE=LOG(EEM1)
      ALFA=1.+ALE/AL10T
      V0=0.18*(4.+ALE/AL10T)*ALFA*(ALFA*VMIN)**0.6666667
      BETA=0.1*(1.+3.*ALE/AL10T)
      B=0.9/(1.+0.4*ALE+0.022*ALE*ALE)
      AA=1.+2.*B*LOG(VC/V0)
      IF(AA.LE.1.) AA=1.05
      A1=1.-AA
      CC=EXP(-0.25*A1*A1/B)
      A1R=1./A1
      C1=VMAX**A1
      C2=VC**A1
C
C     SAMPLE V AND RO
C
  50  CALL GRNDM(RNDM,2)
      R=RNDM(1)
      V=(R*C1+(1.-R)*C2)**A1R
      IF(V.LE.VMIN) GOTO 50
      IF(V.LT.V0) THEN
        SCREJ=CC*((V-VMIN)/(V0-VMIN))**BETA*(V0/V)**A1
      ELSE
        SCREJ=CC*(V0/V)**(A1+B*LOG(V/V0))
      ENDIF
      IF(RNDM(2).GT.SCREJ) GOTO 50
      R0MAX= SCREJ*(1.-EM6/(EEM1**2*(1.-V)))
      CALL GRNDM(RNDM,2)
      R0   = R0MAX*(2.*RNDM(1)-1.)
C
C           Energies
C
      EPP  = V*EEM1
      IF(IPAIR.NE.1)THEN
         NGKINE=0
         DESTEP=DESTEP+EPP
         GO TO 60
      ENDIF
      EPOS = 0.5*EPP*(1.+R0)
      EMIN = EPP-EPOS
C
C           Angles
C
      THETA = AMASS/EEM1
      SINTH = SIN(THETA)
      COSTH = COS(THETA)
      PHI   = TWOPI*RNDM(2)
      COSPHI= COS(PHI)
      SINPHI= SIN(PHI)
C
      CALL GFANG(VECT(4),COSAL,SINAL,COSBT,SINBT,ROTATE)
C
C           Positron
C
      NGKINE = 0
      TPOS   = EPOS-EMASS
      IF(TPOS.GT.CUTELE)THEN
         PPOS  = SQRT((EPOS+EMASS)*TPOS)
         NGKINE= NGKINE+1
         GKIN(1,NGKINE)=PPOS*SINTH*COSPHI
         GKIN(2,NGKINE)=PPOS*SINTH*SINPHI
         GKIN(3,NGKINE)=PPOS*COSTH
         GKIN(4,NGKINE)=EPOS
         GKIN(5,NGKINE)=2.
         TOFD(NGKINE)=0.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
         IF(ROTATE)
     +   CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
      ELSE
         DESTEP=DESTEP+TPOS
         CALL GANNI2
      ENDIF
C
C           Electron
C
      TMIN=EMIN-EMASS
      IF(TMIN.GT.CUTELE)THEN
         PMIN  = SQRT((EMIN+EMASS)*TMIN)
         NGKINE= NGKINE+1
         GKIN(1,NGKINE)=-PMIN*SINTH*COSPHI
         GKIN(2,NGKINE)=-PMIN*SINTH*SINPHI
         GKIN(3,NGKINE)=PMIN*COSTH
         GKIN(4,NGKINE)=EMIN
         GKIN(5,NGKINE)=3.
         TOFD(NGKINE)=0.
         GPOS(1,NGKINE) = VECT(1)
         GPOS(2,NGKINE) = VECT(2)
         GPOS(3,NGKINE) = VECT(3)
         IF(ROTATE)
     +   CALL GDROT(GKIN(1,NGKINE),COSAL,SINAL,COSBT,SINBT)
      ELSE
         DESTEP=DESTEP+TMIN
      ENDIF
C
C           Correct muon track
C
  60  GEKIN  = GEKIN-EPP
      GETOT  = GEKIN+AMASS
      VECT(7)= SQRT((GETOT+AMASS)*GEKIN)
      CALL GEKBIN
C
C           Update probabilities
C
 900  CALL GRNDM(RNDM,1)
      ZINTPA = -LOG(RNDM(1))
      SLPAIR = SLENG
      STEPPA = BIG
C
      END
