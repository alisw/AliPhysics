*
* $Id$
*
* $Log$
* Revision 1.2  1997/11/13 09:25:40  gunter
* Correction By Laszlo Urban; protect against divide by 0. if AL becomes 0.
* around GEKIN 2.93E-5.
* Reported by harald@psiclu.cern.ch (signed by keller@biomed.ee.ethz.ch)
*
* Revision 1.1.1.1  1995/10/24 10:21:22  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/03 02/08/94  19.37.36  by  S.Ravndal
*-- Author :
      SUBROUTINE GBREME
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates discrete hard BREMSSTRAHLUNG by electrons.          *
C.    *                                                                *
C.    *  The secondary photon energy is sampled from a                 *
C.    *  parametrization of the bremsstrahlung calculation of          *
C.    *  Seltzer and Berger (NIM B12,p.95(1985)) for electron          *
C.    *  energies 1 keV - 10 GeV . For higher energies the             *
C.    *  parametrization agrees with the screened Bethe-Heitler        *
C.    *  bremsstrahlung spectrum.Migdal corrections are applied        *
C.    *      by default. The routine works ( together                  *
C.    *  with the routines GBRELE and GBRSGE ) without the Migdal      *
C.    *  corrections using the Patchy switch +USE,BETHE.               *
C.    *                                                                *
C.    *  NOTE :                                                        *
C.    *      BCUTE is the cut-off energy above which the photon energy *
C.    *      spectrum is generated.                                    *
C.    *                                                                *
C.    *    ==>Called by : GTELEC                                       *
C.    *       Authors    L.Urban    29/10/93 *********                 *
C.    *                                                                *
C.    *  13/11/97 bug corrected by L.Urban                             *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcjloc.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcking.inc"
#include "geant321/gcmate.inc"
#include "geant321/gccuts.inc"
      DIMENSION RNDM(2)
      PARAMETER (CORFAC=0.805485E-10)
      PARAMETER (TLIM=0.000999)
      PARAMETER  (
     +AL00=-0.205398E+01,AL01= 0.238815E-01,AL02= 0.525483E-03,
     +AL10=-0.769748E-01,AL11=-0.691499E-01,AL12= 0.222453E-02,
     +AL20= 0.406463E-01,AL21=-0.101281E-01,AL22= 0.340919E-03,
     +A10 = 0.467733E+01,A11 =-0.619012E+00,A12 = 0.202225E-01,
     +A20 =-0.734101E+01,A21 = 0.100462E+01,A22 =-0.320985E-01,
     +A30 = 0.293119E+01,A31 =-0.403761E+00,A32 = 0.125153E-01,
     +BL00= 0.104133E+01,BL01=-0.943291E-02,BL02=-0.454758E-03,
     +BL10= 0.119253E+00,BL11= 0.407467E-01,BL12=-0.130718E-02,
     +BL20=-0.159391E-01,BL21= 0.727752E-02,BL22=-0.194405E-03,
     +B10 = 0.423071E+01,B11 =-0.610995E+00,B12 = 0.195531E-01,
     +B20 =-0.712527E+01,B21 = 0.969160E+00,B22 =-0.274255E-01,
     +B30 = 0.269925E+01,B31 =-0.363283E+00,B32 = 0.955316E-02)
C.    ------------------------------------------------------------------
C.
C             Ensure cut-off avoids infra-red catastrophe.
C
      IF (GEKIN.LE.BCUTE) GO TO 30
      KCASE=NAMEC(9)
C
*******************************>  Z3=Q(JPROB+2)
      Z3=(Z*(Z+1.))**0.3333333
      IF(Z3.LE.0.)GO TO 30
      Z32=Z3**2
      EEL1   = GETOT
      XC=BCUTE/GEKIN
      ALXC=LOG(XC)
      U=LOG(GEKIN/EMASS)
      U2=U**2
      V=LOG(Z)
*
      IF(GEKIN.LE.TLIM) THEN
         AL0=AL00+AL01*Z3+AL02*Z32
         AL1=AL10+AL11*Z3+AL12*Z32
         AL2=AL20+AL21*Z3+AL22*Z32
         AL=AL0+AL1*U+AL2*U2
         BL0=BL00+BL01*Z3+BL02*Z32
         BL1=BL10+BL11*Z3+BL12*Z32
         BL2=BL20+BL21*Z3+BL22*Z32
         BL=BL0+BL1*U+BL2*U2
         GMAX=1.+AL*XC+BL*XC**2
         IF(GEKIN.LT.0.0001) THEN
            G1=1.+AL+BL
            IF(G1.GT.GMAX) GMAX=G1
*
            IF(ABS(AL).GT.1.e-6) THEN
              X0=-BL/(2.*AL)
              IF((XC.LT.X0).AND.(X0.LT.1.)) THEN
                 G0=1.+AL*X0+BL*X0**2
                 IF(G0.GT.GMAX) GMAX=G0
              ENDIF
            ENDIF
*
         ENDIF
      ELSE
         U3=U2*U
         A1=A10+A11*Z3+A12*Z32
         A2=A20+A21*Z3+A22*Z32
         A3=A30+A31*Z3+A32*Z32
         AH=1.+A1/U+A2/U2+A3/U3
         B1=B10+B11*Z3+B12*Z32
         B2=B20+B21*Z3+B22*Z32
         B3=B30+B31*Z3+B32*Z32
         BH=0.75+B1/U+B2/U2+B3/U3
*
         F=4*V-0.55*V**2
         DEL0=136.*EMASS/(Z3*EEL1)
         EPC=XC*GEKIN/EEL1
         DC=DEL0*EPC/(1.-EPC)
         CC=42.392-F
*
         IF(DC.LE.1.) THEN
            DC2=DC**2
            F1=(42.392-7.796*DC+1.961*DC2-F)/CC
            F2=(41.734-6.484*DC+1.250*DC2-F)/CC
         ELSE
            F1=(42.24-8.368*LOG(DC+0.952)-F)/CC
            IF(F1.LT.0.) F1=0.
            F2=F1
         ENDIF
*
         GMAX=(1.-AH*EPC)*F1+BH*EPC**2*F2
      ENDIF
*
      CORR0=CORFAC*DENS*Z/A
      EPM=GEKIN/EEL1
      SC0=1.+CORR0/(EPM*EPM)
*
*     sample photon energy  ( according to 1/Ephoton)
*
   10 CALL GRNDM(RNDM,2)
*
      X=EXP(RNDM(1)*ALXC)
      EP=X*GEKIN/EEL1
*
*     Migdal correction for Ephoton->0. or no correction (Bethe)
*
#if !defined(CERNLIB_BETHE)
      CORR=SC0/(1.+CORR0/(EP*EP))
#endif
#if defined(CERNLIB_BETHE)
      CORR=1.
#endif
*
*     calculate rejection function g(x)
*
      IF(GEKIN.LE.TLIM) THEN
         G=1.+AL*X+BL*X**2
      ELSE
         D=DEL0*EP/(1.-EP)
         IF(D.LE.1.) THEN
            D2=D**2
            F1=(42.392-7.796*D+1.961*D2-F)/CC
            F2=(41.734-6.484*D+1.250*D2-F)/CC
         ELSE
            F1=(42.24-8.368*LOG(D+0.952)-F)/CC
            IF(F1.LT.0.) F1=0.
            F2=F1
         ENDIF
         G=(1.-AH*EP)*F1+BH*EP**2*F2
      ENDIF
      G=G*CORR/GMAX
      IF(RNDM(2).GT.G) GOTO 10
*
*     photon energy is sampled according to the Seltzer-Berger spectrum
*
      EGAMMA=EEL1*EP
C
C        CUT ON ENERGY THRESHOLD ?
C
      IF((IBREM.NE.1).OR.(EGAMMA.LE.CUTGAM)) THEN
         DESTEP = DESTEP + EGAMMA
         GO TO 20
      ENDIF
C
C             Generate emitted photon angles with respect to a Z-axis
C             defined along parent track. PHI is generated isotropically
C             and THETA is assigned a universal angular distribution
C
      EMASS1 = EMASS
      THETA  = GBTETH(EEL1, EMASS1, EP)*EMASS/EEL1
      SINTH  = SIN(THETA)
      COSTH  = COS(THETA)
      CALL GRNDM(RNDM,1)
      PHI    = TWOPI*RNDM(1)
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
C
C             Polar co-ordinates to momentum components.
C
      NGKINE = NGKINE + 1
      GKIN(1,1)=EGAMMA*SINTH*COSPHI
      GKIN(2,1)=EGAMMA*SINTH*SINPHI
      GKIN(3,1)=EGAMMA*COSTH
      GKIN(4,1)=EGAMMA
      GKIN(5,1)=1.
      TOFD(1)  =0.
      GPOS(1,1) = VECT(1)
      GPOS(2,1) = VECT(2)
      GPOS(3,1) = VECT(3)
C
C             Rotate photon into GEANT system
C
      CALL GVROT(VECT(4),GKIN)
C
C             Correct track for lost energy
C
   20 CONTINUE
      GEKIN = GEKIN - EGAMMA
      GETOT = GEKIN + EMASS
      VECT(7)=SQRT (ABS((GETOT+EMASS)*GEKIN))
      CALL GEKBIN
C
C             Update probabilities
C
   30 CALL GRNDM(RNDM,1)
      ZINTBR=-LOG(RNDM(1))
C
      END
