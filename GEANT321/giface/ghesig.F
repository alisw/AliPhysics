*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:15  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.37  by  S.Giani
*-- Author :
      FUNCTION GHESIG(PPART,EKIN,AVER,A,Z,W,KK,DENS,QCOR,LPART)
C
C *** CALCULATION OF THE PROBABILITIES FOR (IN)ELASTIC INTERACTIONS ***
C *** OF STABLE PARTICLES ON PROTON AND NEUTRON                     ***
C *** NVE 07-APR-1988 ***
C
C CALLED BY : GPGHEI
C ORIGIN : F.CARMINATI, H.FESEFELDT (ROUTINE INTACT 06-OCT-1987)
C
C *** IPART DENOTES THE GHEISHA PARTICLE INDEX ***
C
C CONVENTION :
C
C   PARTICLE                 IPART
C   ------------------------------
C   GAMMA                    1
C   NEUTRINO                 2
C   POSITRON                 3
C   ELECTRON                 4
C   MUON +                   5
C   MUON -                   6
C   PION +                   7
C   PION 0                   8
C   PION -                   9
C   KAON +                  10
C   KAON 0 S                11
C   KAON 0 L                12
C   KAON -                  13
C   PROTON                  14
C   PROTON BAR              15
C   NEUTRON                 16
C   NEUTRON BAR             17
C   LAMBDA                  18
C   LAMBDA BAR              19
C   SIGMA +                 20
C   SIGMA 0                 21
C   SIGMA -                 22
C   SIGMA + BAR             23
C   SIGMA 0 BAR             24
C   SIGMA - BAR             25
C   XSI 0                   26
C   XSI -                   27
C   XSI 0 BAR               28
C   XSI - BAR               29
C   DEUTERON                30
C   TRITON                  31
C   ALPHA                   32
C   OMEGA -                 33
C   OMEGA - BAR             34
C   NEW PARTICLES           35
C
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
C
C --- FOR CROSS-SECTION INFORMATION SEE "PCSDATA" ---
C
#include "geant321/gsecti.inc"
C --- GHEISHA COMMONS ---
#include "geant321/s_result.inc"
#include "geant321/s_prntfl.inc"
C
      DIMENSION A(KK),Z(KK),W(KK)
C
      DIMENSION ALPHA(35),ALPHAC(41),IPART2(7),CSA(4)
      DIMENSION PARTEL(35),PARTIN(35),ICORR(35),INTRC(35)
C
      PARAMETER (ONETHR=1./3.)
      SAVE ALPHA,ALPHAC,PARTEL,PARTIN,CSA,IPART2,ICORR,INTRC
C
#include "geant321/s_csdim.inc"
#include "geant321/pcodim.inc"
#include "geant321/s_csdat.inc"
#include "geant321/pcodat.inc"
C
      DATA ALPHA    / 6*0.7,
     +                0.75 ,0.75 ,0.75 ,
     +                0.76,0.76 ,0.76 ,0.76 ,
     +                0.685,0.63 ,0.685,0.63,0.685,0.63,
     +                3*0.685,3*0.63,2*0.685,2*0.63,
     +                3*0.7,0.685,0.63,0.7/
      DATA ALPHAC    /1.2,1.2,1.2,1.15,0.90,0.91,0.98,1.06,1.10,1.11,
     +                1.10,1.08,1.05,1.01,0.985,0.962,0.945,0.932,
     +                0.925,0.920,0.920,0.921,0.922,0.923,0.928,0.931,
     +                0.940,0.945,0.950,0.955,0.958,0.962,0.965,0.976,
     +                0.982,0.988,0.992,1.010,1.020,1.030,1.040/
      DATA PARTEL/6*0.,29*1./
      DATA PARTIN/6*0.,1.00,0.00,1.05,1.20,1.35,1.30,1.20,1.00,1.30,
     +            1.00,1.30,1.00,1.30,1.00,1.00,1.00,1.30,1.30,1.30,
     +            1.00,1.00,1.30,1.30,1.00,1.,1.,1.,1.3,1./
      DATA ICORR /14*1, 0, 1, 0, 1, 0, 3*1, 3*0, 2*1, 2*0, 4*1, 2*0/
      DATA INTRC /6*0, 1, 0, 12*1, 0, 2*2, 0, 10*1, 0/
C
C CROSS SECTIONS ON NUCLEUS ARE KNOWN ONLY FOR PIONS AND PROTONS.
C THE GENERAL LAW SIGMA(A)=1.25*SIGMA(TOT,PROTON)*A**ALPHA IS VALID
C ONLY FOR MOMENTA > 2 GEV.THE PARAMETRIZATION DONE HERE GIVES ONLY
C A BEHAVIOUR AVERAGED OVER MOMENTA AND PARTICLE TYPES.
C FOR A DETECTOR WITH ONLY A FEW MATERIALS IT'S OF COURSE MUCHBETTER
C TO USE TABLES OF THE MEASURED CROSS SECTIONS .
C FOR ELEMENTS WITH THE FOLLOWING ATOMIC NUMBERS MEASURED CROSS
C SECTIONS ARE AVAILABLE (SEE "PCSDATA").
C
C                 H   AL     CU     PB
      DATA  CSA  /1. ,27.00 ,63.54 ,207.19 /
      DATA IPART2/9,8,7,11,10,13,12/
C
C
C --- Initialise GHESIG and switch to GHEISHA particle code ---
      GHESIG=0.0
      IF(LPART.GT.48)GO TO 160
      IPART=KIPART(LPART)
C
C --- No interaction for gammas, neutrinos, electrons, positrons, muons,
C --- neutral pions, neutral sigmas and antisigmas and new particles.
      IF(INTRC(IPART).EQ. 0) GO TO 160
      P=PPART
      EK=EKIN
C
C --- Initialise the cross-sections with 0.0 ---
      DO 10  K=1,KK
         AIEL(K)=0.0
         AIIN(K)=0.0
         AIFI(K)=0.0
         AICA(K)=0.0
   10 CONTINUE
C
      IF ((IPART .GE. 30) .AND. (IPART .LE. 32)) THEN
C
C --- Take geometrical cross sections for inelastic scattering ---
C --- of deuterons, tritons and alphas ---
         IF (IPART .EQ. 30) THEN
            APART=2.0**ONETHR
         ELSEIF (IPART .EQ. 31) THEN
            APART=3.0**ONETHR
         ELSEIF (IPART .EQ. 32) THEN
            APART=4.0**ONETHR
         END IF
         DO 20 K=1,KK
            AIIN(K)=49.0*(APART+A(K)**ONETHR)**2
   20    CONTINUE
         IF (NPRT(9)) THEN
            WRITE(CHMAIL,10000)
            CALL GMAIL(0,0)
         END IF
C
      ELSE IF ((IPART .EQ. 16) .AND. (EK .LE. 0.0327)) THEN
C
C --- Use tables for low energy neutrons ---
C --- get energy bin ---
         JE2=17
         DO 30 J=2,17
            IF (EK .LT. ELAB(J)) THEN
               JE2=J
               GO TO 40
            END IF
   30    CONTINUE
C
   40    JE1=JE2-1
         EKX=MAX(EK,1.0E-9)
         DELAB=ELAB(JE2)-ELAB(JE1)
         DO 70 K=1,KK
C
C --- Get A bin ---
            JA2=15
            DO 50 J=2,15
               IF (A(K) .LT. CNLWAT(J)) THEN
                  JA2=J
                  GO TO 60
               END IF
   50       CONTINUE
C
   60       JA1=JA2-1
            DNLWAT=CNLWAT(JA2)-CNLWAT(JA1)
C
C --- Use linear interpolation or extrapolation by Y=RCE*X+RCA*X+B ---
C
C --- Elastic cross section ---
C --- E interpolation or extrapolation at JA1 ---
            DY=CNLWEL(JA1,JE2)-CNLWEL(JA1,JE1)
            RCE=DY/DELAB
C --- A interpolation or extrapolation at JE1 ---
            DY=CNLWEL(JA2,JE1)-CNLWEL(JA1,JE1)
            RCA=DY/DNLWAT
            B=CNLWEL(JA1,JE1)-RCE*ELAB(JE1)-RCA*CNLWAT(JA1)
            AIEL(K)=RCE*EK+RCA*A(K)+B
C
C --- Inelastic cross section ---
C --- E interpolation or extrapolation at JA1 ---
            DY=CNLWIN(JA1,JE2)-CNLWIN(JA1,JE1)
            RCE=DY/DELAB
C --- A interpolation or extrapolation at JE1 ---
            DY=CNLWIN(JA2,JE1)-CNLWIN(JA1,JE1)
            RCA=DY/DNLWAT
            B=CNLWIN(JA1,JE1)-RCE*ELAB(JE1)-RCA*CNLWAT(JA1)
            AIIN(K)=RCE*EK+RCA*A(K)+B
C
            IZNO=Z(K)+0.01
            AICA(K)=11.12*CSCAP(IZNO)/(EKX*1.0E6)**0.577
   70    CONTINUE
         IF (NPRT(9)) THEN
            WRITE(CHMAIL,10100)
            CALL GMAIL(0,0)
         END IF
      ELSE
C
C --- Use parametrization of cross section data for all other cases ---
C
         IF (NPRT(9)) THEN
            WRITE(CHMAIL,10200)
            CALL GMAIL(0,0)
         END IF
C
C --- Get momentum bin ---
         J=40
         DO 80 I=2,41
            IF (P .LT. PLAB(I)) THEN
               J=I-1
               GO TO 90
            END IF
   80    CONTINUE
C
C --- Start with  cross sections for scattering on free protons ---
C --- use linear interpolation or extrapolation by Y=RC*X+B     ---
   90    DX=PLAB(J+1)-PLAB(J)
C --- Elastic cross section ---
         DY=CSEL(IPART,J+1)-CSEL(IPART,J)
         RC=DY/DX
         B=CSEL(IPART,J)-RC*PLAB(J)
         AIELIN=RC*P+B
C --- Inelastic cross section ---
         DY=CSIN(IPART,J+1)-CSIN(IPART,J)
         RC=DY/DX
         B=CSIN(IPART,J)-RC*PLAB(J)
         AIININ=RC*P+B
         ALPH=ALPHA(IPART)
         IF(IPART.LT.14) THEN
            DY=ALPHAC(J+1)-ALPHAC(J)
            RC=DY/DX
            B=ALPHAC(J)-RC*PLAB(J)
            CORFAC=RC*P+B
            ALPH=ALPH*CORFAC
C
            IPART3=IPART2(IPART-6)
C
C --- Elastic cross section ---
            DY=CSEL(IPART3,J+1)-CSEL(IPART3,J)
            RC=DY/DX
            B=CSEL(IPART3,J)-RC*PLAB(J)
            XSECEL=RC*P+B
C --- Inelastic cross section ---
            DY=CSIN(IPART3,J+1)-CSIN(IPART3,J)
            RC=DY/DX
            B=CSIN(IPART3,J)-RC*PLAB(J)
            XSECIN=RC*P+B
C
         END IF
C
         DO 100 K=1,KK
            AIEL(K)=AIELIN
            AIIN(K)=AIININ
C
            IF (A(K) .GE. 1.5) THEN
C
C --- A-dependence from parametrization ---
               CREL=1.0
               CRIN=1.0
C --- Get medium bin  1=Hydr.  2=Al  3=Cu  4=Pb ---
               I=3
               IF (A(K) .LT. 50.0) I=2
               IF (A(K) .GT. 100.0) I=4
               IF ((IPART .EQ. 14) .OR. (IPART .EQ. 16)) THEN
C
C --- Protons and neutrons ---
C
C --- Elastic cross section ---
                  DY=CSPNEL(I-1,J+1)-CSPNEL(I-1,J)
                  RC=DY/DX
                  B=CSPNEL(I-1,J)-RC*PLAB(J)
                  XSECEL=RC*P+B
C --- Inelastic cross section ---
                  DY=CSPNIN(I-1,J+1)-CSPNIN(I-1,J)
                  RC=DY/DX
                  B=CSPNIN(I-1,J)-RC*PLAB(J)
                  XSECIN=RC*P+B
                  IF (AIEL(K) .GE. 0.001) CREL=XSECEL/(0.36*AIEL(K)*
     +            CSA(I)**1.17)
                  AITOT=AIEL(K)+AIIN(K)
                  IF (AITOT .GE. 0.001) CRIN=XSECIN/(AITOT*CSA(I)**
     +            ALPH)
C
               ELSEIF (IPART .LT. 15) THEN
C
C --- Calculate correction factors from values on Al,Cu,Pb for all ---
C --- mesons use linear interpolation or extrapolation by Y=RC*X+B ---
C --- Note that data is only available for pions and protons
                  WGCH=0.5
                  IF (A(K) .LT. 20.0) WGCH=0.5+0.5*EXP(-(A(K)-1.0))
                  AIEL(K)=WGCH*AIEL(K)+(1.0-WGCH)*XSECEL
                  AIIN(K)=WGCH*AIIN(K)+(1.0-WGCH)*XSECIN
C
C --- This section not for kaons ---
                  IF (IPART .LT. 10) THEN
C
C --- Elastic cross section ---
                     DY=CSPIEL(I-1,J+1)-CSPIEL(I-1,J)
                     RC=DY/DX
                     B=CSPIEL(I-1,J)-RC*PLAB(J)
                     XSPIEL=RC*P+B
C --- Inelastic cross section ---
                     DY=CSPIIN(I-1,J+1)-CSPIIN(I-1,J)
                     RC=DY/DX
                     B=CSPIIN(I-1,J)-RC*PLAB(J)
                     XSPIIN=RC*P+B
C
                     IF (AIEL(K) .GE. 0.001) CREL=XSPIEL/(0.36* AIEL(K)
     +               *CSA(I)**1.17)
                     AITOT=AIEL(K)+AIIN(K)
                     IF (AITOT .GE. 0.001) CRIN=XSPIIN/(AITOT*CSA(I)
     +               **ALPH)
                  END IF
               END IF
               AIIN(K)=CRIN*(AIIN(K)+AIEL(K))*A(K)**ALPH
               AIEL(K)=CREL*0.36*AIEL(K)*A(K)**1.17
               AIEL(K)=AIEL(K)*PARTEL(IPART)
               AIIN(K)=AIIN(K)*PARTIN(IPART)
            END IF
  100    CONTINUE
C
C --- Fission cross sections ---
C --- A-dependence given by  sigma(3 MeV)=-67.0+38.7*Z**(4/3)/A ---
      END IF
      I=21
      DO 110 II=1,21
         IF (EK .LT. EKFISS(II)) THEN
            I=II
            GO TO 120
         END IF
  110 CONTINUE
C
  120 CONTINUE
      DO 130 K=1,KK
C
C --- No fission for materials with A < 230 ---
         IF (A(K) .GE. 230.) THEN
C
C --- Only data for U(233), U(235) and Pu(239) for Ek .le. 0.01 GeV ---
            J=4
            IF (EK .LE. 0.01) THEN
C
C --- Distinguish U(233), U(235), Pu(239) and U(238) and rest by J=1,4
               IF ((Z(K) .EQ. 92.0).AND.(ABS(A(K)-233.0).LT.0.5))THEN
                  J=1
               ELSEIF((Z(K).EQ.92.0).AND.(ABS(A(K)-235.0).LT.0.5))THEN
                  J=2
               ELSEIF((Z(K).EQ.94.0).AND.(ABS(A(K)-239.0).LT.0.5))THEN
                  J=3
               END IF
            END IF
            IF(J.EQ.4) THEN
C
               Z43BA=Z(K)**(4.0*ONETHR)/A(K)
               Z43BA=MAX(-67.0+38.7*Z43BA,0.)
            ELSE
               Z43BA=1.
            ENDIF
C
C --- Energy dependence taken from U(238) ---
C --- approximated as step-function ---
            AIFI(K)=CSFISS(J,I)*Z43BA
         END IF
  130 CONTINUE
C
C --- Corrections for compounds ---
C --- These corrections should only be applied to inorganic scintill. ---
C --- apply the correction only if user selected it within GEANT ---
      IF (QCOR .GT. 0.0.AND.KK.GT.1) THEN
C --- Do not apply corrections for anti-baryons ---
         IF (ICORR(IPART).EQ.1) THEN
C
C --- ACC40 between 0.3 and 0.5 for Pi+ and Pi0, 0.2 for other mesons ---
            ACC40=0.325
C           IF (IPART .GE. 9) ACC40=0.20
C --- ACC40 = 0.15 for baryons ---
            IF (IPART .GE. 14) ACC40=0.15
C --- ACCA=0.08 for all pions, 0.02 for all other particles ---
            ACCA= 0.08
            IF (IPART .GE. 10) ACCA=0.02
C
            ACCB=0.32*(ACC40-ACCA)
            ACC=ACCA-ACCB*LOG(EK)
            IF (ACC .GT. 0.5) ACC=0.5
            IF (ACC .GT. 0.0) THEN
C
               CAVER=AVER**ACC
               CAVER=1.+(CAVER-1.)*QCOR
               DO 140 K=1,KK
                  AIEL(K)=AIEL(K)*CAVER
                  AIIN(K)=AIIN(K)*CAVER
                  AIFI(K)=AIFI(K)*CAVER
                  AICA(K)=AICA(K)*CAVER
  140          CONTINUE
            END IF
         END IF
      END IF
C
C --- Calculate interaction probability ---
C
C --- Correction factor for high (P > 100 GeV/C) energies ---
      CORH=1.0
C --- Assume a LOG(P) dependence of the correction factor with values ---
C P = 100 GeV/C  ==> CORH = 1.
C P =   1 TeV/C  ==> CORH = 1.25
C     DX=LOG(1000.)-LOG(100.)
C     DY=1.25-1.
C     RC=DY/DX
C     B=1.-RC*LOG(100.)
C     CORH=RC*LOG(P)+B
      IF (P .GT. 100.) CORH=0.1085736156*LOG(P)+0.5
      ALAM=0.0
      DO 150 K=1,KK
         AFACT=AVO*1E-3*DENS*W(K)/A(K)
         AIEL(K)=MAX(CORH*AIEL(K)*AFACT,0.)
         AIIN(K)=MAX(CORH*AIIN(K)*AFACT,0.)
         AIFI(K)=MAX(CORH*AIFI(K)*AFACT,0.)
         AICA(K)=MAX(CORH*AICA(K)*AFACT,0.)
C
         ALAM=ALAM+AIEL(K)+AIIN(K)+AIFI(K)+AICA(K)
  150 CONTINUE
C
C --- Pass the interaction probability to GEANT ---
      GHESIG=ALAM
C
      GO TO 999
C
C --- Printout of skipped particles in case of interface debug ---
  160 IF (NPRT(9)) THEN
         WRITE(CHMAIL,10300) IPART
         CALL GMAIL(0,0)
      END IF
10000 FORMAT(' *GHESIG* GEOM X-SECT. FOR INEL. SCAT. OF D,T AND ALPHA')
10100 FORMAT(' *GHESIG* X-SECT. FROM LOW ENERGY NEUTRON TABLES')
10200 FORMAT(' *GHESIG* X-SECT. FROM PARAMETRIZATION OF DATA')
10300 FORMAT('0*GHESIG* GHEISHA PARTICLE ',I3,' SKIPPED')
  999 END
