*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:04  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.44  by  S.Giani
*-- Author :
*$ CREATE NUCREL.FOR
*COPY NUCREL
*
*=== nucrel ===========================================================*
*
C
C   HJM 17/11/88
C
C                MODIFIED ELASTIC SCATTERING ROUTINE
C                (COMP. KMU-HEP INTERNAL NOTE 88-01)
C
C------------------------------------------------------------------
C
      SUBROUTINE NUCREL(IT,PLAB,EKIN,CX,CY,CZ,ANUC,WEE)
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
C********************************************************************
C     VERSION JULY 81 BY             J. RANFT
C                                    LEIPZIG
C     LAST CHANGE A. FERRARI 26-9-89: RECOIL ENERGY SCORING ADDED
C
C     THIS IS A SUBROUTINE OF FLUKA82 TO SAMPLE ELASTIC INTERACTIONS
C
C     INPUT VARIABLES:
C     IT       = TYPE OF THE PRIMARY
C     PLAB     = PRIMARY PARTICLE LAB MOMENTUM (GEV/C)
C     CX,CY,CZ = PRIMARY PARTICLE DIRECTION COSINES
C     ANUC     = ATOMIC WEIGHT OF THE TARGET NUCLEUS
C
C     OUTPUT VARIABLES:
C     CCX,CCY,CCZ = DIRECTION COSINES OF SCATTERED PARTICLE
C
C********************************************************************
C
      PARAMETER (AMUC12 = 0.9314943228D0)
#include "geant321/aadat.inc"
#include "geant321/finuc.inc"
*------ Common part has been added by A. Ferrari ----------------------*
#include "geant321/part2.inc"
C----------------------------------------------------------------
      REAL RNDM(1)
      DIMENSION ITT(NALLWP),SHP(6),BP(6),XNN(6,6),BA(6,6),XNA(6),AM(6)
      SAVE ITT,SHP,BP,XNN,BA,XNA,AM
      DATA ITT/1,2,0,0,0,0,0,1,2,0,0,5,3,4,5,6,1,2,5,1,1,1,3,5,6,
     &         3,0,0,0,1,2,2,2,1,2,1,2,1,2/
C     ITT(IT) - INTERNAL PARTICLE NUMBER
      DATA AM/ 9.01D0, 12.01D0, 26.98D0, 63.55D0, 118.69D0, 207.19D0/
C     AM - ATOMIC WEIGHTS OF SEVERAL TARGET MATERIALS
      DATA XNA/3.5D0, 3.4D0, 4.5D0, 6.7D0, 8.2D0, 9.5D0/
C     XNA(AM) - EFFECTIVE NUCLEON NUMBER
      DATA BP/12.D0, 12.D0, 4*10.D0/
C     BP(IIT) - INCOHERENT SLOPE (GEV/C)**-2
      DATA SHP/2*40.D0, 2*25.D0, 2*20.D0/
C     SHP(IIT) - TOTAL HADRON - NUCLEON CROSS SECTION  (MB)
C(((((((((((((((      BA(IIT,A)    XNN(IIT,A)       (((((((((((((((((
C     BA - COHERENT SLOPE   (GEV/C)**-2
C     XNN - NORMALIZATION FOR PROJECTILE -NUCLEUS SCATTERING  (MB)
      DATA XNN/
     ( 256.6D0,  240.5D0,  181.0D0,  167.1D0,  144.0D0,  139.0D0,
     ( 335.3D0,  355.0D0,  241.0D0,  229.7D0,  181.5D0,  196.0D0,
     ( 755.1D0, 1344.0D0,  506.0D0,  480.8D0,  431.0D0,  403.0D0,
     (1656.0D0, 1886.0D0, 1104.0D0, 1144.0D0,  914.0D0,  939.0D0,
     (3005.0D0, 2559.0D0, 2112.0D0, 2068.0D0, 1723.0D0, 1791.0D0,
     (5123.0D0, 5444.0D0, 4071.0D0, 3391.0D0, 3368.0D0, 3168.0D0/
      DATA BA/
     ( 72.0D0,  69.3D0,  65.2D0,  63.4D0,  60.0D0,  65.6D0,
     ( 72.3D0,  75.9D0,  65.5D0,  63.1D0,  59.5D0,  68.1D0,
     (119.4D0, 128.7D0, 107.3D0, 105.0D0, 105.6D0,  98.8D0,
     (201.0D0, 212.0D0, 188.9D0, 183.0D0, 179.0D0, 178.0D0,
     (311.0D0, 315.0D0, 286.0D0, 277.0D0, 270.0D0, 266.0D0,
     (454.0D0, 448.0D0, 435.0D0, 397.0D0, 416.0D0, 403.0D0/
C
C********************************************************************
C     PARAMETRIZATION
C     DSIGMA/DOMEGA=C*EXP(-DS*THETA**2)+E*EXP(-FS*THETA**2)
C                   COHERENT PART             INCOHERENT PART
C********************************************************************
C
C***
      NP=1
      TV=0.D0
C***
      IF(IT.EQ.30) GO TO 301
      IF(ANUC.LT.0.99D0) THEN
         CXR(1)=CX
         CYR(1)=CY
         CZR(1)=CZ
         TKI(1)=EKIN
         PLR(1)=PLAB
         KPART(1)=IT
         RETURN
      ENDIF
C
C--------------------------------------------------------------------
C         HJM 10/88       ELASTIC SCATTERING INTO 2-BODY FINALSTATE
C                         FOR NUCLEON-PROTON INITIAL STATE
C         AF  9/91        ALSO PBAR,NBAR-PROTON
C                         FINAL STATE PARTICLES IN /FINUC/
      IF( (ANUC.LT.1.5D0) .AND. (IT.LE.2 .OR. IT.EQ.8 .OR. IT.EQ.9)
     &    .AND. (EKIN.LT.3.5D0) ) THEN
         CALL NUPREL(IT,EKIN,PLAB,CX,CY,CZ)
         WEI(1)=WEE
         WEI(2)=WEE
C
         RETURN
      ENDIF
C-------------------------------------------------------------------
C
C********************************************************************
C     1=P,N, 2=AP,AN, 3=PI+, 4=PI-, 5=K+,K0, 6=K-,K0 BAR
C********************************************************************
C
      IIT=ITT(IT)
      IF(IIT.EQ.0)RETURN
      AP=PLAB
      AP2=AP**2
      A=ANUC
      IF(IIT.GE.2) GO TO 101
      IF(AP.GT.20.D0) GO TO 101
C
C********************************************************************
C     FOR PROTONS BELOW 20 GEV/C
C     DSIGMA/DOMEGA=12.5*A**1.63*EXP(-14.5*A**0.66*T)+
C                   17.5*A**0.33*EXP(-10*T)    FOR A<62
C
C     DSIGMA/DOMEGA=50*A**1.33*EXP(-60*A**0.33*T)+
C     20*A**0.4*EXP(-10*T)  FOR A>62
C********************************************************************
C
      A3=A**0.3333333333333333D0
      ATAR=A
      IF (ATAR .GT. 62.D0) GO TO 1
      C=12.5D0*A**1.63D0
      DS=14.5D0*A3*A3*AP2
      E=17.5D0*A3
      GO TO 2
  1   CONTINUE
      C=50.D0*A3*A
      DS=60.D0*A3*AP2
      E=20.D0*A**0.4D0
  2   CONTINUE
      FS=10.D0*AP2
      GO TO 3
C
C********************************************************************
C     FOR PROTONS OVER 20 GEV/C AND OTHER PARTICLES
C********************************************************************
C
  101 CONTINUE
      DO 200 I=1,6
      IF(A.LE.AM(I)) GO TO 201
  200 CONTINUE
      K=6
      GO TO 202
  201 K=I
  202 KK=K+1
      IF(KK.GT.6) KK=K-1
      XNEL=XNN(IIT,K)+(A-AM(K))*(XNN(IIT,KK)-XNN(IIT,K))/(AM(KK)-AM(K))
      C=XNEL**2
      DS=LOG(BA(IIT,KK)/BA(IIT,K))/LOG(AM(KK)/AM(K))
      DS=BA(IIT,K)*(A/AM(K))**DS
      DS=DS*AP2
      E=XNA(K)+(A-AM(K))*(XNA(KK)-XNA(K))/(AM(KK)-AM(K))
      E=E*SHP(IIT)**2
      FS=BP(IIT)*AP2
      GO TO 3
C
C********************************************************************
C     FOR HEAVY IONS:
C
C       FOR APROJ>25 OR (APROJ>10 AND ATARG>100)
C DSIG/DOM = 216 * PTOT**1.86 * (AP**0.7 + AT**0.7)**2.2 *
C          EXP(-16.1 * (AP**1.2 + AT**0.9)**0.8 * T) +
C          0.3 * PTOT**1.86 * (AP**0.7 + AT**0.7)**2.2 * EXP(-23 * T)
C
C       FOR LIGHTER NUCLIDES:
C DSIG/DOM = 78 * PTOT**1.78 * (AP + AT**0.9)**2.1 *
C          EXP(-16.1 * (AP**1.2 + AT**0.9)**0.8 * T) +
C          0.5 * PTOT**1.78 * (AP + AT**0.9)**2.1 * EXP(-30 * T)
C
C       THE PARAMETRIZATION IS OBTAINED BY FITTING TO
C       DATA GENERATED BY THE SOFT-SPHERES MODEL.
C       (REF. CHAUVIN & AL. PHYS.REV.C. 28 1970 (1983))
C       THE PARAMETRIZATION IS QUITE CRUDE AND FAILS AT LARGE
C       ANGLES.
C       THERE IS NO EXPERIMENTAL VERIFICATION FOR THE MODEL AT
C       ENERGIES ABOVE 100 MEV/A
C                                   MIKA HUHTINEN 25.8.88
C********************************************************************
C
  301 CONTINUE
      ITARA=NINT(ANUC)
C     ALWAYS LET THE LIGHTER PARTICLE BE THE PROJECTILE (THEN THE
C     PARAMETRIZATION WILL GIVE BETTER RESULTS)
      IF (ITARA.GE.IPROA) THEN
        PROA=IPROA
        TARA=ITARA
      ELSE
        PROA=ITARA
        TARA=IPROA
      ENDIF
C     TOTAL LAB. MOMENTUM
      PLTOT=PLAB*PROA
C     SQUARE OF CMS MOMENTUM
      AP2=((PLTOT*TARA)/(TARA+PROA))**2
 
      DS=16.1D0*(PROA**1.2D0 + TARA**0.9D0)**0.8D0*AP2
 
      IF (PROA.GT.25.D0) GO TO 302
      IF ((PROA.GT.10.D0).AND.(TARA.GT.100.D0)) GO TO 302
 
      CFF=PLTOT**1.78D0*(PROA + TARA**0.9D0)**2.1D0
      C=78.D0*CFF
      E=0.5D0*CFF
      FS=30.D0*AP2
      GO TO 3
 
  302 CONTINUE
      CFF=PLTOT**1.86D0*(PROA**0.7D0 + TARA**0.7D0)**2.2D0
      C=216.D0*CFF
      E=0.3D0*CFF
      FS=23.D0*AP2
C
C********************************************************************
C     SAMPLE THE ANGLE AND ROTATE TO GET NEW DIRECTION COSINES
C********************************************************************
C
    3 CONTINUE
      EDS=1.D-9
      IF(DS.LT.10.5D0) EDS=EXP(-2.D0*DS)
      EDS2=EDS**2
      DD=C*(1.D0-EDS2)/DS
      EFS=1.D-9
      IF(FS.LT.10.5D0) EFS=EXP(-2.D0*FS)
      EFS2=EFS**2
      FF=E*(1.D0-EFS2)/FS
      DF=2.D0*DS
      CALL GRNDM(RNDM,1)
      VV=RNDM(1)*(DD+FF)
      IF (VV.GT.DD) DF=2.D0*FS
  31  CONTINUE
      CALL GRNDM(RNDM,1)
      V2=LOG(RNDM(1))/DF
      IF (ABS(V2).GT.2.D0) GO TO 31
      V1=V2
      COD=V1+1.D0
      SID=SQRT(-V1*(2.D0+V1))
*  ******* Computing the recoil energy!!!! A. Ferrari, 26-9-89 ******* *
      KP    = IPTOKP (IT)
      APP2  = AAM(KP)**2
      EPROJ = EKIN + AAM(KP)
*  +-------------------------------------------------------------------*
*  |  Hydrogen: generate the recoil proton
      IF ( ANUC .LT. 1.2D+00 ) THEN
         ATT  = AAM (1)
         ATT2 = ATT**2
         ETTT = EPROJ + ATT
         PTTL = PLAB
         PTTT = 0.D+00
         HELP = ATT2 - APP2 * ( 1.D+00 - COD ) * ( 1.D+00 + COD )
*  | Now: it can happen for amproj > amprot on hydrogen.
         IF ( HELP .LT. 0.D+00 ) GO TO 3
         PLAB = PLAB * ( ( APP2  + ATT*EPROJ ) * COD + (EPROJ + ATT)
     &        * SQRT (HELP) ) / ( (EPROJ + ATT)**2 - (PLAB*COD)**2 )
         NP = 2
         TV = 0.D+00
         WEI   (1) = WEE
         WEI   (2) = WEE
         KPART (1) = IT
         KPART (2) = 1
         TKI   (1) = SQRT ( PLAB**2 + APP2 ) - AAM (KP)
         TKI   (2) = EKIN - TKI (1)
         EKIN      = TKI (1)
         PLR   (1) = PLAB
         PLR   (2) = SQRT ( TKI (2) * ( TKI (2) + 2.D+00 * AAM (1) ) )
         CALL SFECFE ( SFE, CFE )
         CALL TTRANS ( CX, CY, CZ, COD, SID, SFE, CFE, CCX, CCY, CCZ )
         CXR (1) = CCX
         CYR (1) = CCY
         CZR (1) = CCZ
         SINTHE  = PLAB / PLR (2) * SID
         COSTHE  = SQRT ( ( 1.D+00 - SINTHE ) * ( 1.D+00 + SINTHE ) )
         CFE     = - CFE
         SFE     = - SFE
         CALL TTRANS ( CX, CY, CZ, COSTHE, SINTHE, SFE, CFE, CXR (2),
     &                 CYR (2), CZR (2) )
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |  "Heavy" nuclei
*  |  the following line to patch el. scatt. of an, ap, lambdas, sigmas
*  |  etc on hydrogen ( Now it should be useless )
      ELSE
         ATT  = MAX ( ANUC * AMUC12, AAM (KP) )
         ATT2 = ATT**2
         HELP  = ATT2 - APP2 * ( 1.D+00 - COD ) * ( 1.D+00 + COD )
         PLAB  = PLAB*( (APP2  + ATT*EPROJ)*COD + (EPROJ + ATT)
     &         * SQRT (HELP) ) / ( (EPROJ + ATT)**2 - (PLAB*COD)**2 )
         TV = EPROJ - SQRT ( PLAB**2 + APP2 )
         EKIN = EKIN - TV
*  | ***** The following lines only for debugging *****
         IF ( PLAB .LE. 0.D+00 ) THEN
            WRITE (LUNOUT,*)' *** Nucrel: plab',PLAB,' ***'
            WRITE (LUNERR,*)' *** Nucrel: plab',PLAB,' ***'
         END IF
*  | ***** End debugging lines *****
      END IF
*  |              Now the recoil energy is put into Tv
*  +-------------------------------------------------------------------*
      CALL SFECFE(SFE,CFE)
      CALL TTRANS(CX,CY,CZ,COD,SID,SFE,CFE,CCX,CCY,CCZ)
C---------------------------------------------------------
C   MODIFICATION HJM 28/10/88
C   NOTE:  HEAVY ION OPTION  NOT CONSISTENTLY TREATED NOW]]]
C
      WEI(1)=WEE
      CXR(1)=CCX
      CYR(1)=CCY
      CZR(1)=CCZ
      TKI(1)=EKIN
      PLR(1)=PLAB
      KPART(1)=IT
      RETURN
      END
