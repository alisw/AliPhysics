*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:03  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 28/03/94  01.30.59  by  S.Giani
*-- Author :
#if defined(CERNLIB_HPUX)
$OPTIMIZE OFF
#endif
*$ CREATE PREPRE.FOR
*COPY PREPRE
*
*=== prepre ===========================================================*
*
      SUBROUTINE PREPRE ( WEE, NNEUT, NPROT, NHOLE, ICYCL )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/eva0.inc"
#include "geant321/fheavy.inc"
#include "geant321/finuc.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      REAL RNDM(2)
      COMMON / FKCOSP / C1ST (3), C2ND (3), LEMISS
      LOGICAL LEMISS
      COMMON / FKCMCY / NPCYCL (20,2), IEVT, LOUT
      COMMON / FKPLOC / IABCOU
*
      DIMENSION ERNCM  (2), EPSMX  (2), AMMAFT (2), AMNAFT (2),
     &          ZAFTR  (2), EEXMNM (2), EEXDLT (2), EEXANW (2),
     &          EPSEEX (2), EPSDLT (2), EPSANY (2), EPSFIX (29),
     &          C (3), UMULEG (0:2), ACOLEG (0:2), NPART (2),
     &          NPRFIX (2)
      LOGICAL LSPREJ, LNUCTS
*
      NPART (1) = NNEUT
      NPART (2) = NPROT
      IF ( LHLFIX ) THEN
         NHLFIX = NHLEXP
         NHOLE  = NHOLE - NHLFIX
         EHLFIX = 0.D+00
         DO 50 IHOLE = 1, NHLFIX
            EHLFIX = EHLFIX + HOLEXP (IHOLE)
   50    CONTINUE
         IF ( NHOLE .GT. 0 ) THEN
            RHOIMP = ( RHOIMP * ( NHOLE + NHLFIX ) - RHOEXP ) / NHOLE
            EKFIMP = ( EKFIMP * ( NHOLE + NHLFIX ) - EKFEXP ) / NHOLE
         ELSE
            RHOIMP = RHOAVE
            EKFIMP = EKFAVE
         END IF
      ELSE
         NHLFIX = 0
         EHLFIX = 0.D+00
      END IF
      NHINI  = NHOLE
      IF ( PTRES .LT. ANGLGB ) THEN
         UMO2  = ERES * ERES
         UMO   = SQRT (UMO2)
         EKRES = 0.D+00
         GAMCM = 1.D+00
         ETAX  = 0.D+00
         ETAY  = 0.D+00
         ETAZ  = 0.D+00
         PXORI = 0.D+00
         PYORI = 0.D+00
         PZORI = 0.D+00
         PCMORI = 0.D+00
         CALL RACO ( CXAXCM, CYAXCM, CZAXCM )
      ELSE
         UMO2  = ( ERES - PTRES ) * ( ERES + PTRES )
         UMO   = SQRT (UMO2)
         EKRES = ERES - UMO
         GAMCM = ERES  / UMO
         ETACM = PTRES / UMO
         ETAX  = PXRES / UMO
         ETAY  = PYRES / UMO
         ETAZ  = PZRES / UMO
         ETAPCM = ETAX * PXORI + ETAY * PYORI + ETAZ * PZORI
         PHELP = EKORI + AM (KPORI) - ETAPCM / ( GAMCM + 1.D+00 )
         PCMSX = PXORI - ETAX * PHELP
         PCMSY = PYORI - ETAY * PHELP
         PCMSZ = PZORI - ETAZ * PHELP
         PCMORI = SQRT ( PCMSX**2 + PCMSY**2 + PCMSZ**2 )
         CXAXCM = PCMSX / PCMORI
         CYAXCM = PCMSY / PCMORI
         CZAXCM = PCMSZ / PCMORI
      END IF
      BBTAR = IBTAR
      ZZTAR = ICHTAR
      TVCMS  = UMO - AMNRES
      EHLFIX = MIN ( TVCMS, EHLFIX )
      TVEFF  = TVCMS - EHLFIX
      IF ( TVEFF .LE. 0.D+00 ) GO TO 7000
      TEMNU = SQRT ( TVCMS / ANOW * ALEVEL  )
      NEXMX = SQRT ( 2.D+00 * ANOW * TVCMS / ALEVEL )
      NPTOT  = NPART (1) + NPART (2)
      IF ( NHOLE + NPTOT + NHLFIX .GE. NEXMX .OR.
     &     NPTOT .GT. NINT (0.5D+00 * ANOW) ) GO TO 7000
      ANPTOT = NPTOT
      ANPROT = NPART (2)
      ANNEUT = NPART (1)
      AVEBIN = ( ( BBTAR - ZZTAR - ACOLL + ZCOLL ) * BNENRG (2)
     &       + ( ZZTAR - ZCOLL ) * BNENRG (1) ) / ( BBTAR - ACOLL )
      IAAFT = IBRES-1
      IZAFT = ICRES
      AAFTR     = IAAFT
      ZAFTR (1) = IZAFT
      AMMAFT(1) = AAFTR * AMUAMU + 1.D-03 * FKENER (AAFTR,ZAFTR(1))
      AMNAFT(1) = AMMAFT (1) - ZAFTR (1) * AMELEC + ELBNDE (IZAFT)
      CALL EEXLVL ( IAAFT, IZAFT, EEXDLT (1), EEXMNM (1), EEXDUM )
      BNDGAV = BNENRG (2)
      EPSMX  (1) = UMO - AMNRES - BNDGAV - EHLFIX
      EPSDLT (1) = UMO - AMNRES - BNDGAV - EEXDLT (1)
      EPSEEX (1) = UMO - AMNRES - BNDGAV - EEXMNM (1)
      IZAFT = ICRES - 1
      ZAFTR (2) = IZAFT
      AMMAFT(2) = AAFTR * AMUAMU + 1.D-03 * FKENER (AAFTR,ZAFTR(2))
      AMNAFT(2) = AMMAFT (2) - ZAFTR (2) * AMELEC + ELBNDE (IZAFT)
      CALL EEXLVL ( IAAFT, IZAFT, EEXDLT (2), EEXMNM (2), EEXDUM )
      BNDGAV = BNENRG (1)
      EPSMX  (2) = UMO - AMNRES - BNDGAV - EHLFIX
      IF ( EPSMX (1) + EPSMX (2) .LT. 2.D+00 * TEMNU ) GO TO 7000
      EPSDLT (2) = UMO - AMNRES - BNDGAV - EEXDLT (2)
      EPSEEX (2) = UMO - AMNRES - BNDGAV - EEXMNM (2)
      IF ( NP .LE. NP0 ) THEN
         IF ( KPORI .EQ. 8 ) THEN
            EEXANW (1) = EEXDLT (1)
            EPSANY (1) = MIN ( EPSDLT (1), EPSMX (1) )
            EEXANW (2) = 0.D+00
            EPSANY (2) = EPSMX  (2)
         ELSE IF ( KPORI .EQ. 1 ) THEN
            EEXANW (1) = 0.D+00
            EPSANY (1) = EPSMX  (1)
            EEXANW (2) = EEXDLT (2)
            EPSANY (2) = MIN ( EPSDLT (2), EPSMX (2) )
         ELSE
            EEXANW (1) = 0.D+00
            EPSANY (1) = EPSMX (1)
            EEXANW (2) = 0.D+00
            EPSANY (2) = EPSMX (2)
         END IF
      ELSE
         EEXANW (1) = 0.D+00
         EPSANY (1) = EPSMX (1)
         EEXANW (2) = 0.D+00
         EPSANY (2) = EPSMX (2)
      END IF
      SIGIN0 = PI * ( R0SIGK * RMASS (IBRES-1) )**2
      NEMISS = 0
      MNUCTS = 0
 1000 CONTINUE
         NPRFIX (1) = 0
         NPRFIX (2) = 0
         JEMFIX = 0
         EPRFIX = 0.D+00
         IF ( .NOT. LEMISS ) ICYCL = ICYCL + 1
         IF ( NPTOT .LE. 0 ) GO TO 4600
         IF ( NNUCTS .GT. MNUCTS ) THEN
            MNUCTS = MNUCTS + 1
            JNUCTS = INUCTS (MNUCTS)
            JEMIS0 = 2 - ABS (KPNUCL(JNUCTS)) / 8
            JEMIS1 = JEMIS0
            JDEMIS = 1
            EKFSAV = EKFIMP
            RHOSAV = RHOIMP
            TMPRHO = 0.001D+00 * RHOAVE
            RHOIMP = MAX ( RHNUCL (JNUCTS), TMPRHO )
            TMPEKF = 0.001D+00 * EKFAVE
            EKFIMP = MAX ( EKFNUC (JNUCTS), TMPEKF )
            LNUCTS = .TRUE.
         ELSE
            LNUCTS = .FALSE.
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. ANNEUT / ANPTOT ) THEN
               JEMIS0 = 1
               JEMIS1 = 2
               JDEMIS = 1
            ELSE
               JEMIS0 = 2
               JEMIS1 = 1
               JDEMIS = -1
            END IF
         END IF
         NEXTOT = NPTOT + NHOLE
         RHORED = ACOLL / BBTAR
         IF ( LGDHPR ) THEN
            EKFCNN = ( ZNOW * EKFCEN (1) + ( ANOW - ZNOW ) * EKFCEN (2)
     &             ) / ANOW
            RHOWEI = 0.5D+00
            IF ( LNUCTS ) THEN
               AEFREF = EKFIMP
               AEFRAV = EKFIMP
               RHONUC = RHOIMP
            ELSE IF ( NPTOT .LE. 2 .AND. ICYCL .EQ. 1 ) THEN
               AHLTOT = MAX ( NHOLE, 1 )
               WEIGH1 = NHINI / AHLTOT
               WEIGH2 = 1.D+00 - WEIGH1
               AEFREF = WEIGH1 * EKFIMP + WEIGH2 * EKFAVE
               RHONUC = WEIGH1 * RHOIMP + WEIGH2 * RHOAVE
               AEFRAV = WEIGH1 * EKFIMP + WEIGH2 * EKFAVE
               RHONUC = RHORED * ( ( AHLTOT + 1.D+00 - RHOWEI )
     &                * RHONUC + RHOWEI * RHOAVE ) / ( AHLTOT + 1.D+00 )
               AEFRAV = ( AHLTOT * AEFRAV + EKFAVE )
     &                / ( AHLTOT + 1.D+00 )
            ELSE
               AHLTOT = MAX ( NHOLE, 1 )
               WEIGH1 = NHINI / AHLTOT
               WEIGH2 = 1.D+00 - WEIGH1
               AEFREF = WEIGH1 * EKFIMP + WEIGH2 * EKFAVE
               RHONUC = WEIGH1 * RHOIMP + WEIGH2 * RHOAVE
               AEFRAV = WEIGH1 * EKFIMP + WEIGH2 * EKFAVE
               RHONUC = RHORED * ( ( AHLTOT + 1.D+00 - RHOWEI )
     &                * RHONUC + RHOWEI * RHOAVE ) / ( AHLTOT + 1.D+00 )
               AEFRAV = ( AHLTOT * AEFRAV + EKFAVE )
     &                / ( AHLTOT + 1.D+00 )
            END IF
            CLAMDI = 1.0D+00
         ELSE
            EKFCNN = AEFRMX
            AEFREF = 0.5D+00 * AEFRMX
            AEFRAV = AEFREF
            RHONUC = 0.5D+00**1.5D+00 * RHONU0
            CLAMDI = 0.5D+00
         END IF
         DO 4500 JEMISS = JEMIS0, JEMIS1, JDEMIS
            IF ( NPART (JEMISS) .LE. 0 ) GO TO 4500
            IEMISS = 3 - JEMISS
            RATEC0 = 2.D+00 * AMNUCL (IEMISS) / PLABRC**3 / PISQ
            BNDGEN = AMNAFT (JEMISS) + AMNUCL (IEMISS) - AMNRES
            IF ( .NOT. LNUCTS .AND. EPRFIX .GT. ANGLGB ) THEN
               TVEFF  = TVEFF - EPRFIX
               NPTOT  = NPTOT - NPRFIX (1) - NPRFIX (2)
               NPART (1) = NPART (1) - NPRFIX (1)
               NPART (2) = NPART (2) - NPRFIX (2)
               JEMFIX = JEMISS
               ERNCM0 = ERNCM (JEMISS)
               EPSMX0 = EPSMX (JEMISS)
               BNDGAV = BNENRG (IEMISS)
               EPSMX (JEMISS) = UMO - AMNRES - BNDGAV - EHLFIX - EPRFIX
            ELSE
               JEMFIX = 0
            END IF
            BNDGAV = BNENRG (IEMISS)
            EPSMIN = BNDGEN - BNDGAV
            BNDHLP = BNENRG (IEMISS)
            IF ( EPSMX (JEMISS) - EPSMIN .LT. TEMNU ) GO TO 4500
            IF ( EHLFIX + EPRFIX .LE. 0.5D+00 * EEXANW (JEMISS) ) THEN
               UUMIN  = TVEFF - 0.5D+00 * ( EPSANY (JEMISS)
     &                + EPSMX (JEMISS) ) - BNDGAV
            ELSE
               UUMIN  = TVEFF - EPSMX (JEMISS) - BNDGAV
            END IF
            UUMAX  = TVEFF - BNDGEN
            BNDUSE = BNENRG (IEMISS)
            IF ( LNUCTS ) THEN
               EKNNUC = ENNUC  (JNUCTS) - AMNUCL (IEMISS) + BNDGEN
     &                - RSTNUC (JNUCTS)
               EPSHLP = BNDHLP - BNDGEN + EPSMIN
               EKNNUC = EKNNUC + EPSHLP
               IF ( EKNNUC .LE. EPSMIN ) THEN
                  GO TO 4500
               ELSE IF ( EKNNUC .LE. 0.D+00 ) THEN
                  DELTAE = 0.5D+00 * ( - EPSMIN - EKNNUC )
                  EKNNUC = EKNNUC + DELTAE
               ELSE
                  DELTAE = 0.D+00
               END IF
               ENNUC (JNUCTS) = EKNNUC + AMNUCL (IEMISS)
               PNUCCO = SQRT ( EKNNUC * ( EKNNUC + 2.D+00
     &                * AMNUCL (IEMISS) ) )
               PHELP = PNUCCO / PNUCL (JNUCTS)
               PCMSX = PXNUCL (JNUCTS) * PHELP
               PCMSY = PYNUCL (JNUCTS) * PHELP
               PCMSZ = PZNUCL (JNUCTS) * PHELP
               ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
               PHELP = ENNUC (JNUCTS) - ETAPCM / ( GAMCM + 1.D+00 )
               PCMSX = PCMSX - ETAX * PHELP
               PCMSY = PCMSY - ETAY * PHELP
               PCMSZ = PCMSZ - ETAZ * PHELP
               PCMST = SQRT ( PCMSX**2 + PCMSY**2 + PCMSZ**2 )
               EPS = GAMCM * ENNUC (JNUCTS) - ETAPCM - AMNUCL (IEMISS)
               IF ( EPS .LE. EPSMIN ) THEN
                  GO TO 4500
               ELSE IF ( EPS .LT. ANGLGB * ENNUC (JNUCTS) ) THEN
                  GO TO 4500
               END IF
               EPS   = EPS - DELTAE
               C (1) = PCMSX / PCMST
               C (2) = PCMSY / PCMST
               C (3) = PCMSZ / PCMST
               EPS = MIN ( EPS, HLFHLF * ( EPSMX (JEMISS)
     &             + EPSANY (JEMISS) ) )
               IF ( JEMISS .EQ. 2 ) THEN
                  FLKCOU = DOST ( 1, ZAFTR (JEMISS) )
                  CCOUL  = DOST ( 3, ZAFTR (JEMISS) )
                  ETHRES = FLKCOU * COULBH * ( ZNOW - 1.D+00 )
     &                   / RMASS (IBRES-1)
               END IF
               FLAMMX = 1.D+00
               NTENT  = 1
               ISAMPL = 0
               LSPREJ = .FALSE.
               GO TO 3000
            END IF
            IF ( UUMIN .GE. UUMAX ) GO TO 4500
            EPSINV = 0.5D+00 * ( UMO2 - ( AMNUCL (IEMISS)
     &             + AMNAFT (JEMISS) )**2 ) / AMNAFT (JEMISS)
            EPSWLL = EPSINV + AEFRAV + BNDUSE + AMNUCL (IEMISS)
            BETWLL = SQRT ( 1.D+00 - AMNUSQ (IEMISS) / EPSWLL**2 )
            EKEWLL = EPSWLL - AMNUCL (IEMISS)
            IF ( JEMISS .EQ. 1 ) THEN
               AALPHA = 0.76D+00 + 2.2D+00 / RMASS (IBRES-1)
               BBETAA = ( 2.22D+00 / RMASS (IBRES-1)**2 - 0.05D+00 )
     &                / AALPHA * 1.D-03
               TMP07  = 0.07D+00
               EPSEFF = MIN ( EPSINV, TMP07 )
               SIGINV = 0.1D+00 * XINNEU ( EPSEFF, ZAFTR (JEMISS),
     &                  BBETAA )
               SFTYMX = 1.1D+00
               GSNGLV = 1.5D+00 * ( ANOW - ZNOW ) / EKFCEN (IEMISS)
               IF ( EKEWLL .LE. 0.04D+00 ) THEN
                  SIGMAP = 3.D+03 * PI / ( 1.206D+03 * EKEWLL + (
     &                   -1.86D+00 + 0.09415D+03 * EKEWLL
     &                   + 0.0001306D+06 * EKEWLL**2 )**2 ) + 1.D+03
     &                   * PI / ( 1.206D+03 * EKEWLL + ( 0.4223D+00
     &                   + 0.13D+03 * EKEWLL )**2 )
                  IF ( EKEWLL .LT. 0.02D+00 ) THEN
                     SIGMAN = 0.3333333333333333D+00 * SIGMAP
                  ELSE
                     SIGMAN = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                      + 42.9D+00
                  END IF
               ELSE
                  SIGMAP = 34.10D+00 / BETWLL**2 - 82.2D+00 / BETWLL
     &                   + 82.2D+00
                  SIGMAN = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                   + 42.9D+00
               END IF
               SIGMAP = 0.1D+00 * SIGMAP
               SIGMAN = 0.1D+00 * SIGMAN
               SGNUNU = ( ( ANOW - ZNOW - ANNEUT ) * SIGMAN
     &                + ( ZNOW - ANPROT ) * SIGMAP )
     &                / ( ANOW - ANPTOT )
            ELSE
               EPSCOU = UMO - AMNAFT (JEMISS) - AMNUCL (IEMISS)
               FLKCOU = DOST ( 1, ZAFTR (JEMISS) )
               CCOUL  = DOST ( 3, ZAFTR (JEMISS) )
               ETHRES = FLKCOU * COULBH * ( ZNOW - 1.D+00 )
     &                / RMASS (IBRES-1)
               TMP07  = 0.07D+00
               EPSEFF = MIN ( EPSINV, TMP07 )
               SIGINV = 1.D-01 * XINPRO ( EPSEFF, ZAFTR (JEMISS),
     &                                    ETHRES )
               SFTYMX = 1.2D+00
               GSNGLV = 1.5D+00 * ZNOW / EKFCEN (IEMISS)
               IF ( EKEWLL .LE. 0.04D+00 ) THEN
                  SIGMAN = 3.D+03 * PI / ( 1.206D+03 * EKEWLL + (
     &                   -1.86D+00 + 0.09415D+03 * EKEWLL
     &                   + 0.0001306D+06 * EKEWLL**2 )**2 ) + 1.D+03
     &                   * PI / ( 1.206D+03 * EKEWLL + ( 0.4223D+00
     &                   + 0.13D+03 * EKEWLL )**2 )
                  IF ( EKEWLL .LT. 0.02D+00 ) THEN
                     SIGMAP = 0.3333333333333333D+00 * SIGMAN
                  ELSE
                     SIGMAP = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                      + 42.9D+00
                  END IF
               ELSE
                  SIGMAN = 34.10D+00 / BETWLL**2 - 82.2D+00 / BETWLL
     &                   + 82.2D+00
                  SIGMAP = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                   + 42.9D+00
               END IF
               SIGMAP = 0.1D+00 * SIGMAP
               SIGMAN = 0.1D+00 * SIGMAN
               SGNUNU = ( ( ANOW - ZNOW - ANNEUT ) * SIGMAN
     &                + ( ZNOW - ANPROT ) * SIGMAP )
     &                / ( ANOW - ANPTOT )
            END IF
            ZITA   = AEFRAV / EKEWLL
            IF ( ZITA .LE. 0.5D+00 ) THEN
               PZITA = 1.D+00 - 1.4D+00 * ZITA
            ELSE
               PZITA = 1.D+00 - 1.4D+00 * ZITA + 0.4D+00 * ZITA
     &               * ( 2.D+00 - 1.D+00 / ZITA )**2.5D+00
            END IF
            RATEC  = SIGINV * RATEC0 * EPSINV / GSNGLV
            IF ( RATEC .LE. 0.D+00 ) GO TO 4500
            ALAMDC = BETWLL / RATEC
            ALAMDI = 1.D+00 / ( SGNUNU * RHONUC * PZITA ) + MAX (TWOTWO
     &             * PI * PLABRC * BETWLL * EPSWLL / AMNUSQ (IEMISS),
     &             RZNUCL )
            AMUTOT = 1.D+00 / ALAMDI + 1.D+00 / ALAMDC
            RATEI  = BETWLL / ALAMDI
            FLAMMX = SFTYMX * RATEC / ( RATEC + CLAMDI * RATEI )
            IF ( FLAMMX .LT. 0.D+00 ) THEN
               GO TO 4500
            END IF
            LSPREJ = .NOT. ( NHOLE .GT. 2 .OR. NPTOT .GT. 4 )
            IF ( LSPREJ ) THEN
               IF ( NPTOT .EQ. 1 ) THEN
                  IF ( NHOLE .EQ. 0 ) THEN
                     EPS = EPSMX (JEMISS)
                     FLAMMX = 1.D+00
                     NTENT  = 1
                     ISAMPL = 0
                     LSPREJ = .FALSE.
                     GO TO 3000
                  ELSE IF ( NHOLE .EQ. 1 ) THEN
                     ISAMPL = 1
                  ELSE IF ( NHOLE .EQ. 2 ) THEN
                     ISAMPL = 2
                  ELSE
                     ISAMPL = 7
                  END IF
               ELSE IF ( NPTOT .EQ. 2 ) THEN
                  IF ( NHOLE .EQ. 0 ) THEN
                     ISAMPL = 7
                  ELSE IF ( NHOLE .EQ. 1 ) THEN
                     ISAMPL = 3
                  ELSE IF ( NHOLE .EQ. 2 ) THEN
                     ISAMPL = 4
                  ELSE
                     ISAMPL = 7
                  END IF
               ELSE IF ( NPTOT .EQ. 3 ) THEN
                  IF ( NHOLE .EQ. 0 ) THEN
                     ISAMPL = 7
                  ELSE IF ( NHOLE .EQ. 1 ) THEN
                     ISAMPL = 5
                  ELSE IF ( NHOLE .EQ. 2 ) THEN
                     ISAMPL = 6
                  ELSE
                     ISAMPL = 7
                  END IF
               ELSE
                  ISAMPL = 7
               END IF
            ELSE
               ISAMPL = 7
            END IF
            GO TO ( 2100, 2200, 2300, 2400, 2500, 2600, 2700 ), ISAMPL
 2100       CONTINUE
               LSPREJ = .FALSE.
               AEFRPA = AEFREF
               UUMAX  = MIN ( AEFRPA, UUMAX )
               AITOT  = ( AEFRPA**3 - ( AEFRPA + UUMIN - UUMAX )**3 )
     &                / AEFRPA**3
               IF ( UUMAX .GT. TVEFF ) THEN
                  AIPRO = MIN ( AITOT, ONEONE ) * FLAMMX
                  AITOT = AITOT * FLAMMX
               ELSE
                  AITOT = AITOT * FLAMMX
                  AIPRO = AITOT
               END IF
            GO TO 2800
 2200       CONTINUE
               LSPREJ = .FALSE.
               AEFRPA = AEFREF
               UUMAX  = MIN ( TWOTWO * AEFRPA, UUMAX )
               IF ( TVEFF .LE. AEFRPA ) THEN
                  EDENOM = TVEFF**2 * 0.5D+00
               ELSE IF ( TVEFF .LE. 2.D+00 * AEFRPA ) THEN
                  EDENOM = ( TVEFF**2 - 2.D+00 * ( TVEFF - AEFRPA )**2 )
     &                   * 0.5D+00
               ELSE
                  EDENOM = AEFRPA**2
               END IF
               UUMN2 = UUMIN**2
               IF ( UUMAX .GT. AEFRPA ) THEN
                  UUDIV = AEFRPA
                  UUDV2 = UUDIV**2
                  AIPR1 = 0.5D+00 * UUDV2
                  AIPR2 = 0.5D+00 * ( UUMAX - UUDIV ) * ( 3.D+00 * UUDIV
     &                  - UUMAX )
                  IF ( UUMIN .LE. UUDIV ) THEN
                     AITO1 = 0.5D+00 * ( UUDV2 - UUMN2 )
                     AITO2 = AIPR2
                     AIHLP = 0.D+00
                  ELSE
                     AITO1 = 0.D+00
                     AIHLP = 0.5D+00 * ( UUMIN - UUDIV ) * ( 3.D+00
     &                     * UUDIV - UUMIN )
                     AITO2 = AIPR2 - AIHLP
                  END IF
               ELSE
                  UUDIV = UUMAX
                  UUDV2 = UUDIV**2
                  AIPR1 = 0.5D+00 * UUDV2
                  AITO1 = 0.5D+00 * ( UUDV2 - UUMN2 )
                  AITO2 = 0.D+00
                  AIPR2 = 0.D+00
               END IF
               AITOT = AITO1 + AITO2
               AITO1 = AITO1 / AITOT
               AITO2 = AITO2 / AITOT
               AITOT = AITOT * FLAMMX / EDENOM
               IF ( UUMAX .GT. TVEFF ) THEN
                  AIPRO = MIN ( ( AIPR1 + AIPR2 ) / EDENOM, ONEONE )
     &                  * FLAMMX
               ELSE
                  AIPRO = ( AIPR1 + AIPR2 ) * FLAMMX / EDENOM
               END IF
            GO TO 2800
 2300       CONTINUE
               LSPREJ = .FALSE.
               AEFRPA = AEFREF
               IF ( TVEFF .LE. AEFRPA ) THEN
                  EDENOM = 0.25D+00 * TVEFF * TVEFF
               ELSE
                  EDENOM = 0.25D+00 * AEFRPA * ( 2.D+00 * TVEFF
     &                   - AEFRPA )
               END IF
               FCHLP = 0.25D+00 * NPART (JEMISS)
               IF ( UUMAX .LE. AEFRPA ) THEN
                  UUDIV = UUMAX
                  AIPR1 = FCHLP * UUDIV**2
                  AITO1 = AIPR1 - FCHLP * UUMIN**2
                  AITO2 = 0.D+00
                  AIPR2 = 0.D+00
               ELSE
                  UUDIV = AEFRPA
                  AIPR1 = FCHLP * UUDIV**2
                  AIPR2 = FCHLP * ( UUMAX - UUDIV ) * UUDIV
                  IF ( UUMIN .GT. UUDIV ) THEN
                     AITO1 = 0.D+00
                     AIHLP = FCHLP * ( UUMIN - UUDIV ) * UUDIV
                     AITO2 = AIPR2 - AIHLP
                  ELSE
                     AIHLP = 0.D+00
                     AITO1 = AIPR1 - FCHLP * UUMIN**2
                     AITO2 = AIPR2
                  END IF
               END IF
               AITOT = AITO1 + AITO2
               AITO1 = AITO1 / AITOT
               AITO2 = AITO2 / AITOT
               IF ( UUMAX .GT. TVEFF ) THEN
                  DDNPAR = NPART(JEMISS)
                  AIPRO = MIN ( ( AIPR1 + AIPR2 ) / EDENOM,
     &                    DDNPAR ) * FLAMMX
               ELSE
                  AIPRO = ( AIPR1 + AIPR2 ) * FLAMMX / EDENOM
               END IF
            GO TO 2800
 2400       CONTINUE
               LSPREJ = .TRUE.
               AEFRPA = AEFREF
               IF ( TVEFF .LE. AEFRPA ) THEN
                  LSPREJ = .FALSE.
                  FEMAX  = 1.D+00
                  ISAMPL = 7
               ELSE IF ( TVEFF .LE. 2.D+00 * AEFRPA ) THEN
                  FEMAX = TVEFF**3 / ( TVEFF**3 - 2.D+00 * ( TVEFF
     &                  - AEFRPA )**3 )
               ELSE
                  FEMAX = 0.1666666666666667D+00 * TVEFF**3
     &                  / ( AEFRPA**2 * ( TVEFF - AEFRPA ) )
               END IF
               UUMNR  = ( UUMIN / TVEFF )**(NEXTOT-1)
               UUMXR  = ( UUMAX / TVEFF )**(NEXTOT-1)
               AITOT  = NPART (JEMISS) * FLAMMX * ( UUMXR - UUMNR )
     &                * FEMAX
               AIPRO  = NPART (JEMISS) * FLAMMX * UUMXR * FEMAX
            GO TO 2800
 2500       CONTINUE
               LSPREJ = .FALSE.
               AEFRPA = AEFREF
               IF ( TVEFF .LE. AEFRPA ) THEN
                  EDENOM = TVEFF**3 / 36.D+00
               ELSE
                  EDENOM = ( TVEFF**3 - ( TVEFF - AEFRPA )**3 )
     &                   / 26.D+00
               END IF
               IF ( UUMAX .GT. AEFRPA ) THEN
                  UUDIV = AEFRPA
                  UUDV3 = UUDIV**3
                  UUMN3 = UUMIN**3
                  FCHLP = NPART(JEMISS) / 36.D+00
                  AIPR1 = FCHLP * UUDV3
                  AIPR2 = 3.D+00 * FCHLP * UUMAX * UUDIV * ( UUMAX
     &                  - UUDIV )
                  IF ( UUMIN .LE. UUDIV ) THEN
                     AITO1 = AIPR1 - FCHLP * UUMN3
                     AITO2 = AIPR2
                     AIHLP = 0.D+00
                  ELSE
                     AITO1 = 0.D+00
                     AIHLP = 3.D+00 * FCHLP * UUMIN * UUDIV * ( UUMIN
     &                     - UUDIV )
                     AITO2 = AIPR2 - AIHLP
                  END IF
               ELSE
                  UUDIV = UUMAX
                  UUDV3 = UUDIV**3
                  UUMN3 = UUMIN**3
                  FCHLP = NPART(JEMISS) / 36.D+00
                  AIPR1 = FCHLP * UUDV3
                  AITO1 = AIPR1 - FCHLP * UUMN3
                  AIPR2 = 0.D+00
                  AITO2 = 0.D+00
               END IF
               AITOT = AITO1 + AITO2
               AITO1 = AITO1 / AITOT
               AITO2 = AITO2 / AITOT
               AITOT = AITOT * FLAMMX / EDENOM
               IF ( UUMAX .GT. TVEFF ) THEN
                  DDNPAR = NPART(JEMISS)
                  AIPRO = MIN ( ( AIPR1 + AIPR2 ) / EDENOM,
     &                    DDNPAR ) * FLAMMX
               ELSE
                  AIPRO = ( AIPR1 + AIPR2 ) * FLAMMX / EDENOM
               END IF
            GO TO 2800
 2600       CONTINUE
               AEFRPA = AEFREF
               IF ( TVEFF .LE. AEFRPA ) THEN
                  EDENOM = TVEFF**4 / 288.D+00
               ELSE IF ( TVEFF .LE. 2.D+00 * AEFRPA ) THEN
                  EDENOM = ( TVEFF**4 - 2.D+00 * ( TVEFF - AEFRPA )**4 )
     &                   / 288.D+00
               ELSE
                  EDENOM = AEFRPA**2 * ( 0.5D+00 * TVEFF**2 - TVEFF
     &                   * AEFRPA + 7.D+00 / 12.D+00 * AEFRPA**2 )
     &                   / 12.D+00
               END IF
               FCHLP = NPART (JEMISS) / 288.D+00
               IF ( UUMAX .GT. 2.D+00 * AEFRPA ) THEN
                  LSPREJ = .TRUE.
                  UUDIV = AEFRPA
                  UUDI2 = 2.D+00 * AEFRPA
                  UUDV4 = UUDIV**4
                  UUMN4 = UUMIN**4
                  UUD24 = UUDI2**4
                  AIPR1 = FCHLP * UUDV4
                  AIPR2 = FCHLP * ( UUD24 - UUDV4 )
                  AIPR3 = 12.D+00 * FCHLP * AEFRPA**2 * UUMAX
     &                  * ( UUMAX - UUDI2 )
                  IF ( UUMIN .LE. UUDIV ) THEN
                     AITO1 = AIPR1 - FCHLP * UUMN4
                     AITO2 = AIPR2
                     AITO3 = AIPR3
                     AIHLP = 0.D+00
                     AIHL2 = 0.D+00
                  ELSE IF ( UUMIN .LE. UUDI2 ) THEN
                     AITO1 = 0.D+00
                     AIHLP = FCHLP * ( UUMN4 - UUDV4 )
                     AITO2 = AIPR2 - AIHLP
                     AITO3 = AIPR3
                     AIHL2 = 0.D+00
                  ELSE
                     AITO1 = 0.D+00
                     AITO2 = 0.D+00
                     AIHL2 = 12.D+00 * FCHLP * AEFRPA**2 * UUMIN
     &                     * ( UUMIN - UUDI2 )
                     AITO3 = AIPR3 - AIHL2
                  END IF
               ELSE IF ( UUMAX .GT. AEFRPA ) THEN
                  LSPREJ = .TRUE.
                  UUDIV = AEFRPA
                  UUDI2 = UUMAX
                  UUDV4 = UUDIV**4
                  UUMN4 = UUMIN**4
                  UUD24 = UUDI2**4
                  AIPR1 = FCHLP * UUDV4
                  AIPR2 = FCHLP * ( UUD24 - UUDV4 )
                  IF ( UUMIN .LE. UUDIV ) THEN
                     AITO1 = AIPR1 - FCHLP * UUMN4
                     AITO2 = AIPR2
                     AIHLP = 0.D+00
                  ELSE
                     AITO1 = 0.D+00
                     AIHLP = FCHLP * ( UUMN4 - UUDV4 )
                     AITO2 = AIPR2 - AIHLP
                  END IF
                  AIPR3 = 0.D+00
                  AITO3 = 0.D+00
               ELSE
                  LSPREJ = .FALSE.
                  UUDI2 = UUMAX
                  UUDIV = UUMAX
                  UUDV4 = UUDIV**4
                  UUMN4 = UUMIN**4
                  AIPR1 = FCHLP * UUDV4
                  AITO1 = AIPR1 - FCHLP * UUMN4
                  AITO2 = 0.D+00
                  AIPR2 = 0.D+00
                  AITO3 = 0.D+00
                  AIPR3 = 0.D+00
               END IF
               AITOT = AITO1 + AITO2 + AITO3
               AITO1 = AITO1 / AITOT
               AITO2 = AITO2 / AITOT
               AITO3 = AITO3 / AITOT
               AITOT = AITOT * FLAMMX / EDENOM
               AIPRO = ( AIPR1 + AIPR2 + AIPR3 ) * FLAMMX / EDENOM
            GO TO 2800
 2700       CONTINUE
               LSPREJ = .FALSE.
               TVLEFF = LOG(TVEFF)
               IF(UUMAX.LE.0.) THEN
                  UULMXR=-100
               ELSE
                  UULMXR = (LOG(UUMAX)-TVLEFF)*(NEXTOT-1)
               ENDIF
               IF(UULMXR.LT.-88) THEN
                  UUMNR = 0.
                  UUMXR = 0.
               ELSE
                  UUMXR = EXP(UULMXR)
                  IF(UUMIN.LE.0.) THEN
                     UULMNR=-100
                  ELSE
                     UULMNR = (LOG(UUMIN)-TVLEFF)*(NEXTOT-1)
                  ENDIF
                  IF(UULMNR.LT.-88) THEN
                     UUMNR = 0.
                  ELSE
                     UUMNR = EXP(UULMNR)
                  ENDIF
               ENDIF
               AITOT  = NPART (JEMISS) * FLAMMX * ( UUMXR - UUMNR )
               AIPRO  = NPART (JEMISS) * FLAMMX * MIN ( UUMXR, ONEONE )
 2800       CONTINUE
            NTENT  = INT  (AIPRO)
            CALL GRNDM(RNDM,1)
            RNTENT = RNDM (1)
            IF ( RNTENT .LT. AIPRO - NTENT ) NTENT = NTENT + 1
 3000       CONTINUE
            ITACC = 0
            DO 4100 IT = 1, NTENT
               GO TO ( 3100, 3200, 3300, 3400, 3500, 3600, 3700 ),
     &            ISAMPL
               GO TO 4000
 3100          CONTINUE
                  CALL GRNDM(RNDM,1)
                  UURND = AEFRPA * ( 1.D+00 - AITOT / FLAMMX
     &                  * RNDM (1) )**0.3333333333333333D+00
                  UURND = AEFRPA - UURND + UUMIN
               GO TO 3800
 3200          CONTINUE
                  CALL GRNDM(RNDM,1)
                  RNDEPS = RNDM (1)
                  IF ( RNDEPS .LT. AITO1 ) THEN
                     RNDEPS = RNDEPS / AITO1
                     UURND  = SQRT ( RNDEPS * ( UUDV2 - UUMN2 ) + UUMN2)
                  ELSE
                     RNDEPS = RNDEPS - AITO1
                     RNDEPS = AIHLP + RNDEPS * AITOT * EDENOM / FLAMMX
                     UURND  = 2.D+00 * UUDIV - SQRT ( UUDIV**2
     &                      - 2.D+00 * RNDEPS )
                  END IF
               GO TO 3800
 3300          CONTINUE
                  CALL GRNDM(RNDM,1)
                  RNDEPS = RNDM (1)
                  IF ( RNDEPS .LT. AITO1 ) THEN
                     RNDEPS = RNDEPS * AITOT
                     UURND  = RNDEPS / FCHLP + UUMIN**2
                     UURND  = SQRT (UURND)
                  ELSE
                     RNDEPS = ( RNDEPS - AITO1 ) * AITOT + AIHLP
                     UURND  = RNDEPS / FCHLP / UUDIV + UUDIV
                  END IF
               GO TO 3800
 3400          CONTINUE
                  CALL GRNDM(RNDM,1)
                  UURND = TVEFF * ( RNDM (1) * AITOT / FLAMMX
     &                  / FEMAX / NPART(JEMISS) + UUMNR )**
     &                  ( 1.D+00 / (NEXTOT-1) )
                  IF ( UURND .LE. AEFRPA ) THEN
                     LSPREJ = .FALSE.
                     FREJE  = 1.D+00
                  ELSE IF ( UURND .LE. 2.D+00 * AEFRPA ) THEN
                     FREJE  = ( UURND**2 - 2.D+00 * ( UURND - AEFRPA )
     &                      **2 ) / UURND**2
                  ELSE
                     FREJE  = 2.D+00 * AEFRPA**2 / UURND**2
                  END IF
               GO TO 3800
 3500          CONTINUE
                  CALL GRNDM(RNDM,1)
                  RNDEPS = RNDM (1)
                  IF ( RNDEPS .LT. AITO1 ) THEN
                     RNDEPS = RNDEPS / AITO1
                     UURND  = RNDEPS * ( UUDV3 - UUMN3 ) + UUMN3
                     UURND  = UURND**0.333333333333333D+00
                     LSPREJ = .FALSE.
                  ELSE IF ( RNDEPS .LT. AITO1 + AITO2 ) THEN
                     RNDEPS = RNDEPS - AITO1
                     RNDEPS = RNDEPS * AITOT * EDENOM / FLAMMX + AIHLP
                     UURND  = 0.5D+00 * ( UUDIV + SQRT ( UUDIV**2
     &                      + 1.333333333333333D+00 * RNDEPS / FCHLP
     &                      / UUDIV ) )
                     LSPREJ = .FALSE.
                  END IF
               GO TO 3800
 3600          CONTINUE
                  CALL GRNDM(RNDM,1)
                  RNDEPS = RNDM (1)
                  IF ( RNDEPS .LT. AITO1 ) THEN
                     RNDEPS = RNDEPS / AITO1
                     UURND  = RNDEPS * ( UUDV4 - UUMN4 ) + UUMN4
                     UURND  = UURND**0.25D+00
                     LSPREJ = .FALSE.
                  ELSE IF ( RNDEPS .LT. AITO1 + AITO2 ) THEN
                     RNDEPS = RNDEPS - AITO1
                     RNDEPS = RNDEPS * AITOT * EDENOM / FLAMMX + AIHLP
                     UURND  = RNDEPS / FCHLP + UUDV4
                     UURND  = UURND**0.25D+00
                     LSPREJ = .TRUE.
                     FREJE  = 1.D+00 - 2.D+00 * ( 1.D+00 - UUDIV / UURND
     &                      )**3
                  ELSE
                     RNDEPS = RNDEPS - AITO1 - AITO2
                     RNDEPS = RNDEPS * AITOT * EDENOM / FLAMMX + AIHL2
                     UURND  = AEFRPA + SQRT ( AEFRPA**2 + RNDEPS /
     &                      ( 12.D+00 * FCHLP * AEFRPA**2 ) )
                     LSPREJ = .FALSE.
                  END IF
               GO TO 3800
 3700          CONTINUE
                  CALL GRNDM(RNDM,1)
                  UURND = TVEFF * ( RNDM (1) * AITOT / FLAMMX
     &                  / NPART(JEMISS) + UUMNR )**
     &                  ( 1.D+00 / (NEXTOT-1) )
 3800          CONTINUE
               IF ( LSPREJ ) THEN
                  CALL GRNDM(RNDM,1)
                  IF ( RNDM (1) .GE. FREJE ) GO TO 4100
               END IF
               EPS    = TVEFF - UURND - BNDGAV
 4000          CONTINUE
               IF ( LPRFIX ) THEN
                  ITACC = ITACC + 1
                  EPSFIX (ITACC) = EPS - EPSMIN + BNDGEN
               END IF
               DEPSEX = EPS - EPSEEX (JEMISS)
               IF ( DEPSEX .GT. 0.D+00 ) THEN
                  IF ( EEXDLT (JEMISS) .LT. EEXMNM (JEMISS) .AND.
     &               DEPSEX .GT. 0.5D+00 * ( EPSDLT (JEMISS) - EPSEEX
     &              (JEMISS) ) ) THEN
                     DEPSEX = EPS - EPSDLT (JEMISS)
                     IF ( DEPSEX .GT. 0.5D+00 * ( EPSDLT (JEMISS)
     &                  - EPSANY (JEMISS) ) ) THEN
                        EPS = EPSANY (JEMISS)
                     ELSE
                        EPS = EPSDLT (JEMISS)
                     END IF
                  ELSE
                     EPS = EPSEEX (JEMISS)
                  END IF
                  IF ( EPS .LE. EPSMIN ) GO TO 4100
               END IF
               AMNHLP = UMO - EPS + EPSMIN - AMNUCL (IEMISS)
               ERNCM (JEMISS) = 0.5D+00 * ( UMO2 + AMNHLP**2
     &                        - AMNUSQ (IEMISS) ) / UMO
               EPS    = UMO - ERNCM (JEMISS) - AMNUCL (IEMISS)
               EPSINV = 0.5D+00 * ( UMO2 - ( AMNUCL (IEMISS)
     &                + AMNHLP )**2 ) / AMNHLP
               EPSWLL = EPSINV + AEFRAV + BNDUSE + AMNUCL (IEMISS)
               BETWLL = SQRT ( 1.D+00 - AMNUSQ (IEMISS) / EPSWLL**2 )
               EKEWLL = EPSWLL - AMNUCL (IEMISS)
               EPSCOU = UMO - AMNHLP - AMNUCL (IEMISS)
               IF ( JEMISS .EQ. 1 ) THEN
                  GSNGLV = 1.5D+00 * ( ANOW - ZNOW ) * SQRT ( ( EKEWLL
     &                   + EKFCEN (IEMISS) - AEFRAV ) / EKFCEN (IEMISS))
     &                   / EKFCEN (IEMISS)
                  AALPHA = 0.76D+00 + 2.2D+00 / RMASS (IBRES-1)
                  BBETAA = ( 2.22D+00 / RMASS (IBRES-1)**2 - 0.05D+00 )
     &                / AALPHA * 1.D-03
                  SIGINV = 0.1D+00 * XINNEU ( EPSINV, ZAFTR (JEMISS),
     &                                        BBETAA )
                  IF ( EKEWLL .LE. 0.04D+00 ) THEN
                     SIGMAP = 3.D+03 * PI / ( 1.206D+03 * EKEWLL + (
     &                      -1.86D+00 + 0.09415D+03 * EKEWLL
     &                      + 0.0001306D+06 * EKEWLL**2 )**2 ) + 1.D+03
     &                      * PI / ( 1.206D+03 * EKEWLL + ( 0.4223D+00
     &                      + 0.13D+03 * EKEWLL )**2 )
                     IF ( EKEWLL .LT. 0.02D+00 ) THEN
                        SIGMAN = 0.3333333333333333D+00 * SIGMAP
                     ELSE
                        SIGMAN = 10.63D+00 / BETWLL**2 - 29.92D+00
     &                         / BETWLL + 42.9D+00
                     END IF
                  ELSE
                     SIGMAP = 34.10D+00 / BETWLL**2 - 82.2D+00 / BETWLL
     &                      + 82.2D+00
                     SIGMAN = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                      + 42.9D+00
                  END IF
                  SIGMAP = 0.1D+00 * SIGMAP
                  SIGMAN = 0.1D+00 * SIGMAN
                  SGNUNU = ( ( ANOW - ZNOW - ANNEUT ) * SIGMAN
     &                   + ( ZNOW - ANPROT ) * SIGMAP )
     &                   / ( ANOW - ANPTOT )
               ELSE
                  GSNGLV = 1.5D+00 * ZNOW * SQRT ( ( EKEWLL
     &                   + EKFCEN (IEMISS) - AEFRAV ) / EKFCEN (IEMISS))
     &                   / EKFCEN (IEMISS)
                  SIGINV = 1.D-01 * XINPRO ( EPSINV, ZAFTR (JEMISS),
     &                                       ETHRES )
                  IF ( EKEWLL .LE. 0.04D+00 ) THEN
                     SIGMAN = 3.D+03 * PI / ( 1.206D+03 * EKEWLL + (
     &                      -1.86D+00 + 0.09415D+03 * EKEWLL
     &                      + 0.0001306D+06 * EKEWLL**2 )**2 ) + 1.D+03
     &                      * PI / ( 1.206D+03 * EKEWLL + ( 0.4223D+00
     &                      + 0.13D+03 * EKEWLL )**2 )
                     IF ( EKEWLL .LT. 0.02D+00 ) THEN
                        SIGMAP = 0.3333333333333333D+00 * SIGMAN
                     ELSE
                        SIGMAP = 10.63D+00 / BETWLL**2 - 29.92D+00
     &                         / BETWLL + 42.9D+00
                     END IF
                  ELSE
                     SIGMAN = 34.10D+00 / BETWLL**2 - 82.2D+00 / BETWLL
     &                      + 82.2D+00
                     SIGMAP = 10.63D+00 / BETWLL**2 - 29.92D+00 / BETWLL
     &                      + 42.9D+00
                  END IF
                  SIGMAP = 0.1D+00 * SIGMAP
                  SIGMAN = 0.1D+00 * SIGMAN
                  SGNUNU = ( ( ANOW - ZNOW - ANNEUT ) * SIGMAN
     &                   + ( ZNOW - ANPROT ) * SIGMAP )
     &                   / ( ANOW - ANPTOT )
               END IF
               ZITA   = AEFRAV / EKEWLL
               IF ( ZITA .LE. 0.5D+00 ) THEN
                  PZITA = 1.D+00 - 1.4D+00 * ZITA
               ELSE
                  PZITA = 1.D+00 - 1.4D+00 * ZITA + 0.4D+00 * ZITA
     &                  * ( 2.D+00 - 1.D+00 / ZITA )**2.5D+00
               END IF
               RATEC  = SIGINV * RATEC0 * EPSINV / GSNGLV
               IF ( RATEC .LE. 0.D+00 ) GO TO 4100
               ALAMDC = BETWLL / RATEC
               ALAMDI = 1.D+00 / ( SGNUNU * RHONUC * PZITA ) + MAX (
     &                  TWOTWO * PI * PLABRC * BETWLL * EPSWLL
     &                / AMNUSQ (IEMISS), RZNUCL )
               AMUTOT = 1.D+00 / ALAMDI + 1.D+00 / ALAMDC
               RATEI  = BETWLL / ALAMDI
               FLAMDA = RATEC / ( RATEC + CLAMDI * RATEI ) / FLAMMX
               CALL GRNDM(RNDM,1)
               IF ( RNDM (1) .LT. FLAMDA ) GO TO 4200
 4100       CONTINUE
            IF ( JEMISS .EQ. JEMIS1 ) GO TO 4500
            IF ( ITACC .LE. NPART (JEMISS) ) THEN
               NPRFIX (JEMISS) = NPRFIX (JEMISS) + ITACC
               DO 4150 IT = 1, ITACC
                  EPRFIX = EPRFIX + EPSFIX (IT)
 4150          CONTINUE
            ELSE
               NPRFIX (JEMISS) = NPART (JEMISS)
               ACCEP = 0.D+00
               DO 4160 IT = 1, ITACC
                  PRACC  = ( NPART (JEMISS) - ACCEP )
     &                   / ( ITACC - IT + 1 )
                  CALL GRNDM(RNDM,1)
                  IF ( RNDM (1) .LT. PRACC ) THEN
                     EPRFIX = EPRFIX + EPSFIX (IT)
                     ACCEP  = ACCEP  + 1.D+00
                  END IF
 4160          CONTINUE
            END IF
            GO TO 4500
 4200       CONTINUE
            ZNOW = ZAFTR (JEMISS)
            ANOW = AAFTR
            IBRES = IBRES - 1
            ICRES = ICRES - JEMISS + 1
            AMMRES = AMMAFT (JEMISS)
            AMNRES = AMNAFT (JEMISS)
            EEPCM  = EPS + AMNUCL (IEMISS)
            ERESCM = UMO - EEPCM
            IF ( LNUCTS ) THEN
            ELSE IF ( ICYCL .EQ. 1 .AND. LEMISS ) THEN
               C (1) = C2ND (1)
               C (2) = C2ND (2)
               C (3) = C2ND (3)
            ELSE IF ( ICYCL .LE. 1 .AND. PCMORI .LE. ANGLGB ) THEN
               C (1) = C1ST (1)
               C (2) = C1ST (2)
               C (3) = C1ST (3)
            ELSE IF ( ICYCL .LE. 1 ) THEN
               PCMSXT = - AMNUCL (IEMISS) * ETAX
               PCMSYT = - AMNUCL (IEMISS) * ETAY
               PCMSZT = - AMNUCL (IEMISS) * ETAZ
               APSHAR = MAX ( ANPTOT + NEMISS - 1, ONEONE )
               PCMSX  = ( CXAXCM * PCMORI + PCMSXT ) / APSHAR - PCMSXT
               PCMSY  = ( CYAXCM * PCMORI + PCMSYT ) / APSHAR - PCMSYT
               PCMSZ  = ( CZAXCM * PCMORI + PCMSZT ) / APSHAR - PCMSZT
*
               PCMSTT = PCMSXT**2 + PCMSYT**2 + PCMSZT**2
               PCMSPT = PCMSX**2  + PCMSY**2  + PCMSZ**2
               ECMSTT = SQRT ( AMNUSQ (IEMISS) + PCMSTT )
               ECMSPT = SQRT ( AMNUSQ (IEMISS) + PCMSPT )
*
               UMMO2  = 2.D+00 * AMNUSQ (IEMISS) + 2.D+00 * ECMSTT
     &                * ECMSPT - 2.D+00 * ( PCMSXT * PCMSX + PCMSYT
     &                * PCMSY + PCMSZT * PCMSZ )
               UMMO   = SQRT ( UMMO2 )
               GAMCCM = ( ECMSTT + ECMSPT ) / UMMO
               ETAXM  = ( PCMSX + PCMSXT ) / UMMO
               ETAYM  = ( PCMSY + PCMSYT ) / UMMO
               ETAZM  = ( PCMSZ + PCMSZT ) / UMMO
               ETAPCM = ETAXM * PCMSX + ETAYM * PCMSY + ETAZM * PCMSZ
               PHELP  = ECMSPT - ETAPCM / ( GAMCCM + 1.D+00 )
               PCCMSX = PCMSX - ETAXM * PHELP
               PCCMSY = PCMSY - ETAYM * PHELP
               PCCMSZ = PCMSZ - ETAZM * PHELP
               PCCMSP = SQRT ( PCCMSX**2 + PCCMSY**2 + PCCMSZ**2 )
               CALL RACO (C(1),C(2),C(3))
               SCAPRO = PCCMSX * C (1) + PCCMSY * C (2) + PCCMSZ * C (3)
               IF ( SCAPRO .LT. 0.D+00 ) THEN
                  C (1) = - C (1)
                  C (2) = - C (2)
                  C (3) = - C (3)
               END IF
               PCCMSX = C(1) * PCCMSP
               PCCMSY = C(2) * PCCMSP
               PCCMSZ = C(3) * PCCMSP
               ETAPCM = ETAXM * PCCMSX + ETAYM * PCCMSY + ETAZM * PCCMSZ
               PHELP  = 0.5D+00 * UMMO + ETAPCM / ( GAMCCM + 1.D+00 )
               PCMSX  = PCCMSX + ETAXM * PHELP
               PCMSY  = PCCMSY + ETAYM * PHELP
               PCMSZ  = PCCMSZ + ETAZM * PHELP
               PCMSPT = SQRT ( PCMSX**2 + PCMSY**2 + PCMSZ**2 )
               C (1)  = PCMSX / PCMSPT
               C (2)  = PCMSY / PCMSPT
               C (3)  = PCMSZ / PCMSPT
               ETAPCM = - ETAPCM
               PHELP  = 0.5D+00 * UMMO + ETAPCM / ( GAMCCM + 1.D+00 )
               PCMSX  = - PCCMSX + ETAXM * PHELP
               PCMSY  = - PCCMSY + ETAYM * PHELP
               PCMSZ  = - PCCMSZ + ETAZM * PHELP
               PCMSPT = SQRT ( PCMSX**2 + PCMSY**2 + PCMSZ**2 )
               C2ND (1) = PCMSX / PCMSPT
               C2ND (2) = PCMSY / PCMSPT
               C2ND (3) = PCMSZ / PCMSPT
            ELSE
               UMULEG (0) = 1.0D+00
               UMULEG (1) = 0.6666666666666667D+00
               UMULEG (2) = 0.25D+00
               ACOLEG (0) = 0.5D+00
               ACOLEG (1) = 1.5D+00 * UMULEG (1)**ICYCL
               ACOLEG (2) = 2.5D+00 * UMULEG (2)**ICYCL
               COSTHE = COSLEG ( ACOLEG )
               SINTHE = SQRT ( ( 1.D+00 - COSTHE ) * ( 1.D+00 + COSTHE )
     &                )
 4300          CONTINUE
                  CALL GRNDM(RNDM,2)
                  RPHI1 = 2.D+00 * RNDM (1) - 1.D+00
                  RPHI2 = 2.D+00 * RNDM (2) - 1.D+00
                  RPHI12 = RPHI1 * RPHI1
                  RPHI22 = RPHI2 * RPHI2
                  RSQ = RPHI12 + RPHI22
               IF ( RSQ .GT. 1.D+00 ) GO TO 4300
               SINPHI = 2.D+00 * RPHI1 * RPHI2 / RSQ
               COSPHI = ( RPHI12 - RPHI22 ) / RSQ
               SINT02 = ( 1.D+00 - CZAXCM ) * ( 1.D+00 + CZAXCM )
               IF ( SINT02 .LT. ANGLSQ ) THEN
                  C (1) = COSPHI * SINTHE
                  C (2) = SINPHI * SINTHE
                  C (3) = CZAXCM * COSTHE
               ELSE
                  SINTH0 = SQRT ( SINT02 )
                  UPRIME = SINTHE * COSPHI
                  VPRIME = SINTHE * SINPHI
                  COSPH0 = CXAXCM / SINTH0
                  SINPH0 = CYAXCM / SINTH0
                  C (1) = UPRIME * COSPH0 * CZAXCM - VPRIME * SINPH0
     &                  + COSTHE * CXAXCM
                  C (2) = UPRIME * SINPH0 * CZAXCM + VPRIME * COSPH0
     &                  + COSTHE * CYAXCM
                  C (3) = - UPRIME * SINTH0 + COSTHE * CZAXCM
               END IF
            END IF
            PCMS  = SQRT ( ( EEPCM - AMNUCL (IEMISS) ) * ( EEPCM
     &            + AMNUCL (IEMISS) ) )
            UMO2   = ( ERESCM - PCMS ) * ( ERESCM + PCMS )
            UMO    = SQRT (UMO2)
            PCMSX = PCMS * C (1)
            PCMSY = PCMS * C (2)
            PCMSZ = PCMS * C (3)
            ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
            NP = NP + 1
            KPART (NP) = 1 - 7 * ( JEMISS - 2 )
            TKI   (NP) = GAMCM * EEPCM + ETAPCM - AMNUCL (IEMISS)
            PHELP = ETAPCM / ( GAMCM + 1.D+00 ) + EEPCM
            WEI   (NP) = WEE
            PLBPX = PCMSX + ETAX * PHELP
            PLBPY = PCMSY + ETAY * PHELP
            PLBPZ = PCMSZ + ETAZ * PHELP
            PHELP = SQRT ( PLBPX * PLBPX + PLBPY * PLBPY
     &            + PLBPZ * PLBPZ )
            CXR   (NP) = PLBPX / PHELP
            CYR   (NP) = PLBPY / PHELP
            CZR   (NP) = PLBPZ / PHELP
            PLR   (NP) = PHELP
            IGREYP = IGREYP + JEMISS - 1
            IGREYN = IGREYN + 2 - JEMISS
            IF ( JEMISS .EQ. 1 ) THEN
               ISTORE = IGREYN
            ELSE
               ISTORE = IGREYP
            END IF
            IF ( LNUCTS ) THEN
               NPCYCL (ISTORE,JEMISS) = -ICYCL
            ELSE
               NPCYCL (ISTORE,JEMISS) = ICYCL
            END IF
            EINTR  = EINTR  + TKI (NP) + AMNUCL (IEMISS)
            PXINTR = PXINTR + CXR (NP) * PLR (NP)
            PYINTR = PYINTR + CYR (NP) * PLR (NP)
            PZINTR = PZINTR + CZR (NP) * PLR (NP)
            IBINTR = IBINTR + IBAR ( KPART (NP) )
            ICINTR = ICINTR + ICH  ( KPART (NP) )
            ERES  = GAMCM * ERESCM - ETAPCM
            EKRES = ERES - UMO
            PHELP = - ETAPCM / ( GAMCM + 1.D+00 ) + ERESCM
            PXRES = - PCMSX + ETAX * PHELP
            PYRES = - PCMSY + ETAY * PHELP
            PZRES = - PCMSZ + ETAZ * PHELP
            PTRES2= PXRES * PXRES + PYRES * PYRES + PZRES * PZRES
            PTRES = SQRT (PTRES2)
            TVCMS = UMO - AMNRES
            IF ( TVCMS .LT. 0.D+00 ) THEN
               IF ( TVCMS .LT. - GAMMIN ) THEN
                  WRITE (LUNOUT,*)' PREPRE: TVCMS',TVCMS
                  WRITE (LUNERR,*)' PREPRE: TVCMS',TVCMS
                  WRITE (77,*)' ^^^ PREPRE: TVCMS',TVCMS
               END IF
               TVCMS = 0.D+00
            END IF
            EHLFIX = MIN ( TVCMS, EHLFIX )
            TVEFF = TVCMS - EHLFIX
            IF ( TVEFF .LE. 0.D+00 ) GO TO 7000
            GAMCM = ERES  / UMO
            ETACM = PTRES / UMO
            ETAX  = PXRES / UMO
            ETAY  = PYRES / UMO
            ETAZ  = PZRES / UMO
            ETAPCM = ETAX * PXORI + ETAY * PYORI + ETAZ * PZORI
            PHELP = EKORI + AM (KPORI) - ETAPCM / ( GAMCM + 1.D+00 )
            PCMSX = PXORI - ETAX * PHELP
            PCMSY = PYORI - ETAY * PHELP
            PCMSZ = PZORI - ETAZ * PHELP
            PCMORI = PCMSX**2 + PCMSY**2 + PCMSZ**2
            IF ( PCMORI .LE. ANGLGB ) THEN
               PCMORI = 0.D+00
               CALL RACO ( CXAXCM, CYAXCM, CZAXCM )
            ELSE
               PCMORI = SQRT ( PCMORI )
               CXAXCM = PCMSX / PCMORI
               CYAXCM = PCMSY / PCMORI
               CZAXCM = PCMSZ / PCMORI
            END IF
            IAAFT = IBRES-1
            IZAFT = ICRES
            AAFTR = IAAFT
            ZAFTR (1) = IZAFT
            AMMAFT(1) = AAFTR * AMUAMU + 1.D-03
     &                * FKENER (AAFTR,ZAFTR(1))
            AMNAFT(1) = AMMAFT (1) - ZAFTR (1) * AMELEC
     &                + ELBNDE (IZAFT)
            CALL EEXLVL ( IAAFT, IZAFT, EEXDLT (1), EEXMNM (1), EEXDUM )
            BNDGAV = BNENRG (2)
            EPSMX  (1) = UMO - AMNRES - BNDGAV - EHLFIX
            EPSDLT (1) = UMO - AMNRES - BNDGAV - EEXDLT (1)
            EPSEEX (1) = UMO - AMNRES - BNDGAV - EEXMNM (1)
            IZAFT = ICRES - 1
            ZAFTR (2) = IZAFT
            AMMAFT(2) = AAFTR * AMUAMU + 1.D-03
     &                * FKENER (AAFTR,ZAFTR(2))
            AMNAFT(2) = AMMAFT (2) - ZAFTR (2) * AMELEC
     &                + ELBNDE (IZAFT)
            CALL EEXLVL ( IAAFT, IZAFT, EEXDLT (2), EEXMNM (2), EEXDUM )
            BNDGAV = BNENRG (1)
            EPSMX  (2) = UMO - AMNRES - BNDGAV - EHLFIX
            EPSDLT (2) = UMO - AMNRES - BNDGAV - EEXDLT (2)
            EPSEEX (2) = UMO - AMNRES - BNDGAV - EEXMNM (2)
            EEXANW (1) = 0.D+00
            EPSANY (1) = EPSMX (1)
            EEXANW (2) = 0.D+00
            EPSANY (2) = EPSMX (2)
            NPART (JEMISS) = NPART (JEMISS) - 1
            IF ( JEMFIX .GT. 0 ) THEN
               NPART (1) = NPART (1) + NPRFIX (1)
               NPART (2) = NPART (2) + NPRFIX (2)
            END IF
            NEMISS = NEMISS + 1
            LEMISS = .TRUE.
            TEMNU = SQRT ( TVCMS / ANOW * ALEVEL  )
            IF ( EPSMX (1) + EPSMX (2) .LT. 2.D+00 * TEMNU )
     &         GO TO 7000
            SIGIN0 = PI * ( R0SIGK * RMASS (IBRES-1) )**2
            NEXMX = SQRT ( 2.D+00 * ANOW * TVCMS / ALEVEL )
            ANPROT =  NPART (2)
            ANNEUT =  NPART (1)
            NPTOT  = NPART (1) + NPART (2)
            ANPTOT = NPTOT
            GO TO 6000
 4500    CONTINUE
 4600    CONTINUE
         IF ( LNUCTS ) THEN
            LEMISS = .TRUE.
            GO TO 6000
         ELSE
            LEMISS = .FALSE.
            NPART (1) = NPART (1) + NPRFIX (1)
            NPART (2) = NPART (2) + NPRFIX (2)
            IF ( JEMFIX .GT. 0 ) THEN
               TVEFF = TVCMS - EHLFIX
               EPSMX (JEMFIX) = EPSMX0
               ERNCM (JEMFIX) = ERNCM0
            END IF
         END IF
 5000    CONTINUE
         ANPROT =  NPART (2)
         ANNEUT =  NPART (1)
         IF ( NPTOT .GT. 0 ) THEN
            PNPROT = ( ZNOW - ANPROT ) * ( 3.D+00 * ANNEUT + ANPROT )
     &             / ( ANPROT * ( ZNOW - ANPROT + 3.D+00 * ( ANOW
     &             - ANNEUT - ZNOW ) ) + ANNEUT * ( 3.D+00 * ( ZNOW
     &             - ANPROT ) + ANOW - ANNEUT - ZNOW ) )
         ELSE
            PNPROT = ZNOW / ANOW
         END IF
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. PNPROT ) THEN
            NPART (2) = NPART (2) + 1
            ZCOLL = ZCOLL - 1.D+00
         ELSE
            NPART (1) = NPART (1) + 1
         END IF
         NPTOT  = NPART (1) + NPART (2)
         ANPTOT = NPTOT
         NHOLE  = NHOLE + 1
         ACOLL  = ACOLL - 1.D+00
         AVEBIN = ( ( BBTAR - ZZTAR - ACOLL + ZCOLL ) * BNENRG (2)
     &          + ( ZZTAR - ZCOLL ) * BNENRG (1) ) / ( BBTAR - ACOLL )
 6000    CONTINUE
         IF ( LNUCTS ) THEN
            RHOIMP = RHOSAV
            EKFIMP = EKFSAV
         END IF
      IF ( NHOLE + NHLFIX + NPTOT .LT. NEXMX .AND.
     &     NPTOT .LE. NINT (0.5D+00 * ANOW) ) GO TO 1000
 7000 CONTINUE
*=== End of subroutine prepre =========================================*
      RETURN
      END
#if defined(CERNLIB_HPUX)
$OPTIMIZE ON
#endif
