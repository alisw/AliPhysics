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
      SUBROUTINE GMUNU
C
C *** GENERATION OF MUON-NUCLEUS INTERACTIONS ***
C *** NVE 16-MAR-1988 CERN GENEVA ***
C
C CALLED BY : GTMUON
C ORIGIN : H.FESEFELDT (ROUTINE CASMU)
C
C ***Revised Sep-90 by C. CHIERA, E. LAMANNA:
C ***Rebinning of vectors TETAL and XEML
C ***to avoid big angle for outgoing muons
C
C  Rebinning reinstated as original by H-J Trost. The correction of
C  the angle of the outgoing muon should take care of the anomalies
C  at large angles. 19-JUN-91.
C
C.    *  This routine is a straigth translation of the Gheisha routine *
C.    *  CASMU in the Geant dialect. The inelastic cross section is    *
C.    *  taken as 0.0003 mb for E<30 GeV, and is slowly increasing for *
C.    *  E>30 GeV. The energy and angle of final state muon is         *
C.    *  generated according to the 'free quark parton model'. Instead *
C.    *  of the virtual photon a real pion is written on working       *
C.    *  common in order to make use of the cascading routines for     *
C.    *  hadron production.                                            *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gsecti.inc"
#include "geant321/gcmulo.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
C
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION TOTPRO
#endif
      DIMENSION VMUOUT(3),COEF(200),XEML(23),TETAL(35)
      DIMENSION RNDM(3)
      LOGICAL FIRST
      SAVE COEF, FIRST, TETAL, XEML
      DATA FIRST /.TRUE./
      DATA TETAL /
     A   1.0000000,  0.9999995,  0.9999990,  0.9999981,  0.9999962,
     A   0.9999943,  0.9999905,  0.9999847,  0.9999752,  0.9999599,
     A   0.9999352,  0.9998951,  0.9998302,  0.9997253,  0.9995556,
     A   0.9992810,  0.9988368,  0.9981183,  0.9969561,  0.9950773,
     A   0.9920409,  0.9871377,  0.9792297,  0.9665010,  0.9460785,
     A   0.9134827,  0.8618938,  0.7813507,  0.6583430,  0.4770452,
     A   0.2247237, -0.0955139, -0.4461272, -0.7495149, -0.9900000 /
      DATA XEML /1.,.998,.997,.996,.995,.994,.992,.99,.97,.95,
     +   .92,.89,.85,.8,.75,.7,.6,.5,.4,.3,.2,.1,.05/
C             COEF contains the value of the integral of the
C             double differential cross section ds/d(e1) d(cost)
C             for the production of the outgoing muon. These
C             values are obtained using the function
C             GMUSIG and are used to normalize the random value
C             used to sample the energy and angle of the outgoing
C             muon.
C
      IF(FIRST) THEN
*
*   Integrate the double differential cross section
*
         DO 8 ICOEF=1, NEKBIN+1
            EINIT=ELOW(ICOEF)+EMMU
            TOTPRO=0.0
            DO 5 ICOST=2,35
               COST = (TETAL(ICOST)+TETAL(ICOST-1))*0.5
               DO 3 IEFIN=2,23
                  EFINAL = (XEML(IEFIN)+XEML(IEFIN-1))*0.5*EINIT
                  TOTPRO = TOTPRO + GMUSIG(EINIT,EFINAL,COST)*
     +                     (TETAL(ICOST-1)-TETAL(ICOST))*
     +                     (XEML(IEFIN-1)-XEML(IEFIN))*EINIT
   3           CONTINUE
   5        CONTINUE
            COEF(ICOEF) = TOTPRO
   8     CONTINUE
         FIRST=.FALSE.
      ENDIF
C
      KCASE  = NAMEC(21)
      IF(VECT(7).LT.0.01) GO TO 9999
C
C               Generate 4-momentum of final state muon
C
C --- IW2TRY loop introduced to avoid W2<0. (HJT/NVE 27-sep-1990) ---
C --- In case not successfull within 100 tries ==> No change made ---
      IW2TRY=0
 10   CONTINUE
      IF (IW2TRY .GT. 100) GO TO 9999
      IW2TRY=IW2TRY+1
C
      TOTPRO = 0.0
      CALL GRNDM(RNDM,1)
      RAN    = RNDM(1)
      HMAX   = (1.-GEKRAT)*COEF(IEKBIN)+GEKRAT*COEF(IEKBIN+1)
      DO 14 ICOST=2,35
        COST   = (TETAL(ICOST)+TETAL(ICOST-1))*.5
        DO 13 IE1=2,23
          E1     = (XEML(IE1)+XEML(IE1-1))*.5*GETOT
          TOTPRO = TOTPRO+GMUSIG(GETOT,E1,COST)
     *       *(TETAL(ICOST-1)-TETAL(ICOST))
     *       *(XEML(IE1-1)-XEML(IE1))*GETOT
          IF(RAN*HMAX.LT.TOTPRO) GOTO 15
  13    CONTINUE
        IE1    = 23
  14  CONTINUE
      ICOST= 35
  15  CALL GRNDM(RNDM,3)
      TETA = ACOS(TETAL(ICOST-1))+
     *     RNDM(1)*(ACOS(TETAL(ICOST))-ACOS(TETAL(ICOST-1)))
      COST = COS(TETA)
      E1   = (XEML(IE1)+RNDM(2)*(XEML(IE1-1)-XEML(IE1)))*GETOT
      IF(E1.LT.EMMU) E1=EMMU+0.0001
      P1   = SQRT(ABS((E1-EMMU)*(E1+EMMU)))
      IF (ABS(COST) .GT. 1.) COST=SIGN(1.,COST)
C
C --- Check value of W2 and in case negative ==> try again ---
      W2=PMASS*(PMASS+2.*(GETOT-E1))-
     +    2.*(GETOT*E1-VECT(7)*P1*COST-EMMU**2)
      IF (W2 .LT. 0.) GO TO 10
C
      SINT   = SQRT(ABS(1.-COST*COST))
      PHI    = TWOPI*RNDM(3)
      VMUOUT(1) = P1*SINT*COS(PHI)
      VMUOUT(2) = P1*SINT*SIN(PHI)
      VMUOUT(3) = P1*COST
C
C               Store muon on stack and write virtual photon on
C                     result common, rotate muon momenta
C
      CALL GVROT(VECT(4),VMUOUT)
      IF(IMUNU.EQ.1) THEN
C
C            Now compute momentum of the outgoing pion
C
        VECT(4) =  VECT(4)*VECT(7)-VMUOUT(1)
        VECT(5) =  VECT(5)*VECT(7)-VMUOUT(2)
        VECT(6) =  VECT(6)*VECT(7)-VMUOUT(3)
        VECT(7) =  SQRT(VECT(4)*VECT(4)+VECT(5)*VECT(5)+VECT(6)*VECT(6))
        VECT(4) =  VECT(4)/VECT(7)
        VECT(5) =  VECT(5)/VECT(7)
        VECT(6) =  VECT(6)/VECT(7)
C
C            Select pi+ or pi-
C
        IPOLD   = IPART
        CALL GRNDM(RNDM,1)
        IPART   = 8.5+RNDM(1)
        JPA     = LQ(JPART-IPART)
        DO 16 J=1,5
  16    NAPART(J) = IQ(JPA+J)
        ITRTYP    = Q(JPA+6)
        AMASS     = Q(JPA+7)
        CHARGE    = Q(JPA+8)
        TLIFE     = Q(JPA+9)
        GETOT     = SQRT(VECT(7)*VECT(7)+AMASS*AMASS)
        GEKIN     = GETOT-AMASS
        CALL GEKBIN
C
C             Now force interaction of the pion
C
C             This piece of code useful only if using the
C             Gheisha-Geant interface
C
        STEPHA = BIG
        IHOLD     = IHADR
        IF(IHADR.NE.3) IHADR     = 1
#if !defined(CERNLIB_USRJMP)
        CALL GUPHAD
#endif
#if defined(CERNLIB_USRJMP)
        CALL JUMPT0(JUPHAD)
#endif
        KK = Q(JMA+11)
        DO 17 K=1,KK
C
C             Forbid elastic scattering
C
          ALAM    = ALAM - AIEL(K)
          AIEL(K) = 0.0
  17    CONTINUE
        NMOLD     = NMEC
#if !defined(CERNLIB_USRJMP)
        CALL GUHADR
#endif
#if defined(CERNLIB_USRJMP)
        CALL JUMPT0(JUHADR)
#endif
        IHADR     = IHOLD
        NMEC      = NMOLD
        ISTOP     = 0
        IPART     = IPOLD
        JPA       = LQ(JPART-IPART)
        DO 26 J=1,5
  26    NAPART(J) = IQ(JPA+J)
        ITRTYP    = Q(JPA+6)
        AMASS     = Q(JPA+7)
        CHARGE    = Q(JPA+8)
        TLIFE     = Q(JPA+9)
      ELSE
        DESTEP = DESTEP+GETOT-SQRT(P1*P1+AMASS*AMASS)
      ENDIF
C
C            Now just put the muon back on the current stack
C
      VECT(7)   =
     + SQRT(VMUOUT(1)*VMUOUT(1)+VMUOUT(2)*VMUOUT(2)+VMUOUT(3)*VMUOUT(3))
      VECT(4)   = VMUOUT(1)/VECT(7)
      VECT(5)   = VMUOUT(2)/VECT(7)
      VECT(6)   = VMUOUT(3)/VECT(7)
      GETOT     = SQRT(VECT(7)*VECT(7)+AMASS*AMASS)
      GEKIN     = GETOT-AMASS
      CALL GEKBIN
C
C             Update probabilities
C
  90  SLMUNU=SLENG
      ZINTMU=GARNDM(DUM)
      STEPMU=BIG
C
 9999 CONTINUE
C
      RETURN
      END
