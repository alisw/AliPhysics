*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 02/07/94  17.22.04  by  S.Giani
*-- Author :
      SUBROUTINE FEREVV ( NHAD, KP, KT, PM, EKM, TXX, TYY, TZZ, ATEMP,
     &                    ZTEMP, IVVFLG )
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  Ferevt90 by A. Ferrari                                              *
*                                                                      *
*      Last change  on  15-aug-93   by   Alfredo Ferrari, INFN - Milan *
*                                                                      *
*      Ferevt calculates hadron nucleon collisions                     *
*      including the Fermi momentum of the target:                     *
*      now there are two entries, one for valence collisions (or for   *
*      sea collisions with one of the last two nucleons) and one       *
*      (fersea) for sea collisions.                                    *
*                                                                      *
*      Nhad = number of final hadrons                                  *
*      Kp and Kt are indices of the projectile and target nucleons     *
*      Pm   = momentum of the projectile (for ferevv entry)            *
*      Ekm  = kinetic energy of the projectile (for both entries)      *
*      Txx, Tyy, Tzz = direction cosines of the incident projectile,   *
*                      THEY ARE CHANGED in the routine for sea intera- *
*                      ctions                                          *
*      Atemp, Ztemp  = mass and atomic number of the residula nucleus  *
*                      after the "use" of the Kt nucleon               *
*      Kprim = index of the real projectile (only for Fersea)          *
*      Eprim = energy of the real projectile after the emission of     *
*              the virtual meson (only for Fersea)                     *
*      Pprim = momentum of the real projectile after the emission of   *
*              the virtual meson (only for Fersea)                     *
*      Eprold= energy of the real projectile before the emission of    *
*              the virtual meson (only for Fersea)                     *
*      Pprold= momentum of the real projectile before the emission of  *
*              the virtual meson (only for Fersea)                     *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/depnuc.inc"
#include "geant321/hadpar.inc"
#include "geant321/inpdat2.inc"
#include "geant321/nucdat.inc"
#include "geant321/parevt.inc"
#include "geant321/part2.inc"
#include "geant321/resnuc.inc"
      COMMON /FKPRIN/ IPRI, INIT
      LOGICAL LSEACL, LSMPAN
      DIMENSION PTHRSH (NALLWP), IJNUCR (NALLWP)
      COMMON /FKEVNT/ LNUCRI, LHADRI
      LOGICAL LNUCRI, LHADRI
      REAL FRNDM(3)
      REAL FRUUN(1)
      SAVE PTHRSH, IJNUCR
      DATA PTHRSH / 16*5.D+00,2*2.5D+00,5.D+00,3*2.5D+00,8*5.D+00,
     &              9*2.5D+00 /
      DATA IJNUCR / 16*1,2*0,1,3*0,8*1,9*0 /
      IVFLAG = IVVFLG
      AMPROJ = AAM (KP)
      AMTAR  = AAM (KT)
      PSPENT = PM
      ESPENT = EKM + AMPROJ
      PSPTOT = PSPENT
      ESPTOT = ESPENT
      LSEACL = .FALSE.
      IF ( AMPROJ .LT. 1.D+00 ) THEN
         AFACT = 2.8D+00
      ELSE
         AFACT = 1.8D+00
      END IF
      URMIN2 = ( AMPROJ + AMTAR )**2 + ( 1.2D+00 + DBLE (IVFLAG)
     &       * AFACT ) * ETHSEA * AMTAR
      GO TO 100
         ENTRY FERSEA ( NHAD, KP, KT, EKM, TXX, TYY, TZZ, ATEMP, ZTEMP,
     &                  KPRIM, EPRIM, PPRIM, EPROLD, PPROLD )
   50    CONTINUE
         LSEACL = .TRUE.
         IVFLAG = 0
         AMPROJ = AAM (KP)
         AMTAR  = AAM (KT)
         PSPTOT = PPROLD
         ESPTOT = EPROLD
         PSPENT = PPROLD - PPRIM
         ESPENT = EKM + AMPROJ
         IF ( KP .EQ. 26 .AND. KT .EQ. 8 ) THEN
            URMIN2 = ( AAM (13) + AAM (56) )**2 + 0.2D+00 * ETHSEA
     &             * AAM (56)
         ELSE IF ( KP .EQ. 23 .AND. KT .EQ. 1 ) THEN
            URMIN2 = ( AAM (14) + AAM (53) )**2 + 0.2D+00 * ETHSEA
     &             * AAM (53)
         ELSE
            URMIN2 = ( AMPROJ + AMTAR )**2 + 1.2D+00 * ETHSEA * AMTAR
         END IF
 100  CONTINUE
      ITJ = MIN ( 2, KT )
      B1SAVE = B1BAMJ
      B2SAVE = B2BAMJ
      IF ( LLASTN ) THEN
         IF ( LLAST1 ) THEN
            PXNUC = PXXINC
            PYNUC = PYYINC
            PZNUC = PZZINC
            EKNUC = EKINC
         ELSE
            LLAST1 = .TRUE.
            EKNUC  = EKLAST
            PXNUC  = PXLAST
            PYNUC  = PYLAST
            PZNUC  = PZLAST
         END IF
         EM    = EKM   + AMPROJ
         EFER  = EKNUC + AMTAR
         EKFER = EFER  - AMNUCL (ITJ) + EBNDNG (ITJ)
         EFRM  = EFRM  + EKFER
         PXFRM = PXFRM + PXNUC
         PYFRM = PYFRM + PYNUC
         PZFRM = PZFRM + PZNUC
         ECHCK  = EM + EFER
         PXCHCK = PXNUC + PM * TXX
         PYCHCK = PYNUC + PM * TYY
         PZCHCK = PZNUC + PM * TZZ
         UMO  = SQRT ( ECHCK**2 - PXCHCK**2 - PYCHCK**2 - PZCHCK**2 )
      ELSE
         CALL GRNDM(FRNDM,3)
         P2  = MAX ( FRNDM (1), FRNDM (2), FRNDM (3) )
         P2  = PFRMMX (ITJ) * P2
         P2SQ   = P2 * P2
         AMTEMP = ATEMP * AMUC12
         AMTMSQ = AMTEMP * AMTEMP
         EKRECL = 0.5D+00 * P2SQ / AMTEMP * ( 1.D+00 - 0.25D+00 * P2SQ
     &          / AMTMSQ )
         EKREC0 = EKRECL
         CALL POLI   (POLC, POLS)
         CALL SFECFE (SFE , CFE )
         CXTA = POLS * CFE
         CYTA = POLS * SFE
         CZTA = POLC
         PXNUC = P2 * CXTA
         PYNUC = P2 * CYTA
         PZNUC = P2 * CZTA
         PXFRM = PXFRM + PXNUC
         PYFRM = PYFRM + PYNUC
         PZFRM = PZFRM + PZNUC
         EFER  = SQRT ( AMNUSQ (ITJ) + PXNUC**2 + PYNUC**2 + PZNUC**2 )
         FRSURV = ANOW /  DBLE ( IBTAR )
         IF ( .NOT. LSEACL .OR. FRSURV .LT. 0.66D+00 .OR. ATEMP .LT.
     &        40.D+00 ) THEN
            IATEMP = NINT (ATEMP)
            IZTEMP = NINT (ZTEMP)
            AMMRES = AMUAMU * ATEMP + 1.D-03 * FKENER ( ATEMP, ZTEMP )
            AMMRE2 = AMMRES * AMMRES
            EKRECL = SQRT ( P2SQ + AMMRE2 ) - AMMRES
            EKREC0 = EKRECL
            ELEFT  = ETTOT  - EINTR  - EUZ - ESPTOT - EFER + EKRECL
     &             + V0WELL (ITJ) + EBNDNG (ITJ)
            PXLAST = PXTTOT - PXINTR - PUX - PSPTOT * TXX - PXNUC
            PYLAST = PYTTOT - PYINTR - PUY - PSPTOT * TYY - PYNUC
            PZLAST = PZTTOT - PZINTR - PUZ - PSPTOT * TZZ - PZNUC
            PPLAS2 = PXLAST**2 + PYLAST**2 + PZLAST**2
            IUMO = 0
            DELEFT = AMMTAR - AMNTAR - ( DBLE (ICHTAR) - ZTEMP )
     &             * AMELEC
 200        CONTINUE
               EELEFT = ELEFT + DELEFT
               UMO2 = EELEFT**2 - PPLAS2
               DELTU2 = AMMRE2 - UMO2
               IF ( DELTU2 .GT. 0.D+00 ) THEN
                  IUMO = IUMO + 1
                  IF ( LSEACL ) THEN
                     DELTAE = 0.5D+00 * DELTU2 / EELEFT
                     IF ( IUMO .GT. 2 ) THEN
                     ELSE IF ( DELTAE .LT. 2.D+00 * EKREC0 ) THEN
                        EKRECL = EKRECL + DELTAE
                        ELEFT  = ELEFT  + DELTAE
                        GO TO 200
                     ELSE
                        EKRECL = EKRECL + EKREC0
                        ELEFT  = ELEFT  + EKREC0
                        EELEFT = ELEFT + DELEFT
                        UMO2 = EELEFT**2 - PPLAS2
                        DELTU2 = AMMRE2 - UMO2
                        PFDTPL = PXNUC * PXLAST + PYNUC * PYLAST
     &                         + PZNUC * PZLAST
                        PPLAST = SQRT ( PPLAS2 )
                        DELTPR = 0.51D+00 * DELTU2 / ( PPLAST - EELEFT
     &                         * PFDTPL / ( EFER * PPLAST ) )
                        DELTPR = SIGN ( MIN ( ABS ( DELTPR ), HLFHLF
     &                         * P2 ), DELTPR ) / PPLAST
                        PXXINC = PXLAST * DELTPR
                        PYYINC = PYLAST * DELTPR
                        PZZINC = PZLAST * DELTPR
                        PXLAST = PXLAST - PXXINC
                        PYLAST = PYLAST - PYYINC
                        PZLAST = PZLAST - PZZINC
                        PXFRM = PXFRM + PXXINC
                        PYFRM = PYFRM + PYYINC
                        PZFRM = PZFRM + PZZINC
                        PXNUC = PXNUC + PXXINC
                        PYNUC = PYNUC + PYYINC
                        PZNUC = PZNUC + PZZINC
                        EFER0 = EFER
                        EFER  = SQRT ( PXNUC**2 + PYNUC**2 +
     &                                 PZNUC**2 + AMNUSQ  (ITJ) )
                        DELTAE = EFER0 - EFER
                        IF ( DELTAE .GT. 0.D+00 ) THEN
                           ELEFT = ELEFT + DELTAE
                        ELSE
                           EFER = EFER0
                        END IF
                        PPLAS2 = PXLAST**2 + PYLAST**2 + PZLAST**2
                        DELTAE = EKREC0
                        GO TO 200
                     END IF
                  ELSE
                     DELTAE = 0.5D+00 * DELTU2 / EELEFT
                     IF ( IUMO .GT. 3 ) THEN
                        WRITE ( LUNOUT, * )' Ferevv: valence call,',
     &                  ' impossible to get',
     &                  ' enough invariant mass for the residual',
     &                  ' nucleus!!',IATEMP,IZTEMP,AMMRE2,UMO2,UMO,
     &                    DELTAE
                        WRITE ( LUNERR, * )' Ferevv: valence call,',
     &                  ' impossible to get',
     &                  ' enough invariant mass for the residual',
     &                  ' nucleus!!',IATEMP,IZTEMP,AMMRE2,UMO2,UMO,
     &                    DELTAE
                     ELSE IF ( DELTAE .LT. 3.D+00 * EKREC0 ) THEN
                        EKRECL = EKRECL + DELTAE
                        ELEFT  = ELEFT  + DELTAE
                        GO TO 200
                     ELSE
                        EKRECL = EKRECL + 2.D+00 * EKREC0
                        ELEFT  = ELEFT  + 2.D+00 * EKREC0
                        EELEFT = ELEFT + DELEFT
                        UMO2 = EELEFT**2 - PPLAS2
                        DELTU2 = AMMRE2 - UMO2
                        PPDTPL = PM * ( PXLAST * TXX + PYLAST * TYY
     &                         + PZLAST * TZZ )
                        PPLAST = SQRT ( PPLAS2 )
                        DELTPR = 0.51D+00 * DELTU2 / ( PPLAST - EELEFT
     &                         * PPDTPL / ( ( EKM + AMPROJ ) * PPLAST
     &                           ) )
                        TMPPM  = 0.3D+00 * PM
                        TMPPPL = 0.8D+00 * PPLAST
                        DELTPR = SIGN ( MIN ( ABS ( DELTPR ), TMPPM
     &                          , TMPPPL ), DELTPR )
     &                         / PPLAST
                        PXXINC = PXLAST * DELTPR
                        PYYINC = PYLAST * DELTPR
                        PZZINC = PZLAST * DELTPR
                        PXLAST = PXLAST - PXXINC
                        PYLAST = PYLAST - PYYINC
                        PZLAST = PZLAST - PZZINC
                        PXFRM = PXFRM + PXXINC
                        PYFRM = PYFRM + PYYINC
                        PZFRM = PZFRM + PZZINC
                        PXXINC = PM * TXX + PXXINC
                        PYYINC = PM * TYY + PYYINC
                        PZZINC = PM * TZZ + PZZINC
                        PM  = SQRT ( PXXINC**2 + PYYINC**2
     &                             + PZZINC**2 )
                        TXX = PXXINC / PM
                        TYY = PYYINC / PM
                        TZZ = PZZINC / PM
                        DELTAE = EKM
                        EKM = SQRT ( PM * PM + AMPROJ * AMPROJ )
     &                      - AMPROJ
                        DELTAE = DELTAE - EKM
                        ELEFT  = ELEFT  + DELTAE
                        EFRM   = EFRM   - DELTAE
                        PPLAS2 = PXLAST**2 + PYLAST**2 + PZLAST**2
                        PSPENT = PM
                        ESPENT = EKM + AMPROJ
                        PSPTOT = PSPENT
                        ESPTOT = ESPENT
                        GO TO 200
                     END IF
                  END IF
                  UMO = ( ESPENT + EFER - V0WELL (ITJ) - EKRECL
     &                  - EBNDNG (ITJ) )**2
     &                - ( TXX * PSPENT + PXNUC )**2
     &                - ( TYY * PSPENT + PYNUC )**2
     &                - ( TZZ * PSPENT + PZNUC )**2
                  IF ( UMO .LT. URMIN2 ) THEN
                     WRITE ( LUNOUT, * )' Ferevv: impossible to get',
     &                      ' enough invariant mass for interaction',
     &                        IATEMP,IZTEMP,AMMRE2,UMO2,UMO
                     WRITE ( LUNERR, * )' Ferevv: impossible to get',
     &                      ' enough invariant mass for interaction',
     &                        IATEMP,IZTEMP,AMMRE2,UMO2,UMO
                  END IF
               END IF
            IF ( IUMO .GT. 0 .AND. .NOT. LSEACL .AND. DELTU2 .LT.
     &           0.D+00 ) THEN
               DELTAE = SQRT ( DELTU2 + EELEFT**2 )
     &                - EELEFT
               DELTAE = DELTAE + 10.D+00 * ANGLGB * AMMRES
               EKRECL = EKRECL + DELTAE
               ELEFT  = ELEFT  + DELTAE
               EELEFT = ELEFT + DELEFT
               UMO2 = EELEFT**2 - PPLAS2
               DELTU2 = AMMRE2 - UMO2
            END IF
         ELSE
            AMTEMP = ATEMP * AMUC12
            AMTMSQ = AMTEMP * AMTEMP
            EKRECL = 0.5D+00 * P2SQ / AMTEMP * ( 1.D+00 - 0.25D+00
     &             * P2SQ / AMTMSQ )
            EKREC0 = EKRECL
         END IF
         EKFER = EFER - AMNUCL (ITJ)
         IF ( LSEACL ) THEN
            AMPRI2 = EPRIM * EPRIM - PPRIM * PPRIM
            PXCHCK = PXNUC + PPROLD * TXX
            PYCHCK = PYNUC + PPROLD * TYY
            PZCHCK = PZNUC + PPROLD * TZZ
            P2CHCK = PXCHCK**2 + PYCHCK**2 + PZCHCK**2
            ECHCK  = EPROLD + EFER
            UMO2   = ECHCK**2 - P2CHCK
            IF ( UMO2 .LT. 0.D+00 ) THEN
               LRESMP = .TRUE.
               WRITE (LUNERR,*)' FEREVV: SEA INT. UMO2 < 0 ',UMO2
               RETURN
            END IF
            UMO    = SQRT ( UMO2 )
            RMIN2  = URMIN2 / UMO2
            GAMCM = ECHCK  / UMO
            ETAX  = PXCHCK / UMO
            ETAY  = PYCHCK / UMO
            ETAZ  = PZCHCK / UMO
            ETACM = SQRT ( ETAX**2 + ETAY**2 + ETAZ**2 )
            CXCMS = ETAX / ETACM
            CYCMS = ETAY / ETACM
            CZCMS = ETAZ / ETACM
            CALL SFECFE ( SFE, CFE )
            IF ( ABS (CZCMS) .GT. 1.D-04 ) THEN
               CXXTR  = - SFE * CZCMS
               CYYTR  = CFE * CZCMS
               CZZTR  = CXCMS * SFE - CYCMS * CFE
            ELSE IF ( ABS (CYCMS) .GT. 1.D-04 ) THEN
               CXXTR  = CYCMS * CFE
               CYYTR  = CZCMS * SFE - CXCMS * CFE
               CZZTR  = - SFE * CYCMS
            ELSE
               CXXTR  = CYCMS * SFE - CZCMS * CFE
               CYYTR  = - SFE * CXCMS
               CZZTR  = CFE * CXCMS
            END IF
            TXXOLD = TXX
            TYYOLD = TYY
            TZZOLD = TZZ
            PCMSMX = PPRIM - ETACM * ( EPRIM - ETACM * PPRIM / ( GAMCM
     &             + 1.D+00 ) )
            RMAX2  = 1.D+00 + AMPRI2 / UMO2 - 2.D+00 * ( EPRIM -
     &               ETACM * PCMSMX ) / ECHCK
            IF ( RMAX2 .GT. RMIN2 ) THEN
               PTRANS = 0.D+00
            ELSE
               PTRANS = 0.D+00
            END IF
 
1220        CONTINUE
               PLONG2 = ( PPRIM - PTRANS ) * ( PPRIM + PTRANS )
               PLONGI = SQRT ( PLONG2 )
               PXLAST = PTRANS * CXXTR + PLONGI * CXCMS
               PYLAST = PTRANS * CYYTR + PLONGI * CYCMS
               PZLAST = PTRANS * CZZTR + PLONGI * CZCMS
               PCMSLN = PLONGI - ETACM * ( EPRIM - ETACM * PLONGI /
     &                ( GAMCM + 1.D+00 ) )
               RURM2  = 1.D+00 + AMPRI2 / UMO2 - 2.D+00 * ( EPRIM -
     &                  ETACM * PCMSLN ) / ECHCK
               TXX = PXLAST / PPRIM
               TYY = PYLAST / PPRIM
               TZZ = PZLAST / PPRIM
               IF ( RURM2 .LE. RMIN2 ) THEN
                  PTRANS = 0.5D+00 * PTRANS
                  WRITE ( LUNERR, * )' Ferevv: R2 < Rmin2 for Pt',
     &                                 RURM2, RMIN2, URMIN2
                  IF ( PTRANS .GT. ANGLGB ) GO TO 1220
               END IF
            PCMSX  = PCMSLN * CXCMS + PTRANS * CXXTR
            PCMSY  = PCMSLN * CYCMS + PTRANS * CYYTR
            PCMSZ  = PCMSLN * CZCMS + PTRANS * CZZTR
            ERCMS  = 0.5D+00 * ( UMO2 * ( 1.D+00 + RURM2 ) - AMPRI2 )
     &             / UMO
            ETAPCM = PCMSLN * ETACM
            ECHCK  = GAMCM * ERCMS - ETAPCM
            PHELP  = - ETAPCM / (GAMCM + 1.D+00) + ERCMS
            PXCHCK = - PCMSX + ETAX * PHELP
            PYCHCK = - PCMSY + ETAY * PHELP
            PZCHCK = - PCMSZ + ETAZ * PHELP
            ECHCK = ECHCK - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
            UMO   = ECHCK**2 - PXCHCK**2 - PYCHCK**2
     &            - PZCHCK**2
            IF ( UMO .LT. 0.D+00 ) THEN
               LRESMP = .TRUE.
               RETURN
            END IF
            UMO   = SQRT ( UMO )
         ELSE
            EM = EKM + AMPROJ - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
            ECHCK  = EM + EFER
            PXCHCK = PXNUC + PM * TXX
            PYCHCK = PYNUC + PM * TYY
            PZCHCK = PZNUC + PM * TZZ
            UMO    = ECHCK**2 - PXCHCK**2 - PYCHCK**2
     &             - PZCHCK**2
            IF ( UMO .LT. 0.D+00 ) THEN
               WRITE (LUNOUT,*)' *** Ferevv: Umo < 0 ',UMO
               WRITE (LUNERR,*)' *** Ferevv: Umo < 0 ',UMO
               LRESMP = .TRUE.
               NHAD = 0
               RETURN
            END IF
            UMO    = SQRT ( UMO )
         END IF
         EFRM  = EFRM  + EKFER - V0WELL (ITJ) - EKRECL
         TVEUZ = TVEUZ + V0WELL (ITJ) - EKFER + EKRECL
      END IF
      UMO2   = UMO * UMO
      EPROJX = HLFHLF * ( UMO2 - AMPROJ**2 - AMTAR**2 ) / AMTAR
      IF ( EPROJX .LE. AMPROJ ) THEN
         WRITE (LUNOUT,*)' Ferevv: Eprojx < Amproj after kin. sel. !',
     &                     EPROJX, AMPROJ, LSEACL
         WRITE (LUNERR,*)' Ferevv: Eprojx < Amproj after kin. sel. !',
     &                     EPROJX, AMPROJ, LSEACL
         PXFRM = PXFRM - PXNUC
         PYFRM = PYFRM - PYNUC
         PZFRM = PZFRM - PZNUC
         EFRM  = EFRM  - EKFER + V0WELL (ITJ) + EKRECL
         TVEUZ = TVEUZ - V0WELL (ITJ) + EKFER - EKRECL
         IF ( LSEACL ) GO TO 50
         LRESMP = .TRUE.
         NHAD   = 0
         RETURN
      END IF
      PPROJX = SQRT ( ( EPROJX - AMPROJ ) * ( EPROJX + AMPROJ ) )
      ETOTX  = EPROJX + AMTAR
      PTOSCA = PXCHCK * TXX + PYCHCK * TYY + PZCHCK * TZZ
      PXTART = PXCHCK - PTOSCA * TXX
      PYTART = PYCHCK - PTOSCA * TYY
      PZTART = PZCHCK - PTOSCA * TZZ
      PTRASQ = PXTART**2 + PYTART**2 + PZTART**2
      AMTRSQ = AMTAR**2  + PTRASQ
      UMOTR2 = UMO2 + PTRASQ
      UMOTR  = SQRT (UMOTR2)
      PPARSQ = ECHCK**2 - UMOTR2
      PPARTT = SQRT (PPARSQ)
      GAMCMS = ECHCK  / UMOTR
      ETACMS = PPARTT / UMOTR
      EPRCMS = HLFHLF * ( UMOTR2 + AMPROJ**2 - AMTRSQ ) / UMOTR
      PPRCMS = SQRT ( ( EPRCMS - AMPROJ ) * ( EPRCMS + AMPROJ ) )
      EPRLAB = GAMCMS * EPRCMS + ETACMS * PPRCMS
      ETRLAB = ECHCK  - EPRLAB
      PPRLAB = SQRT ( ( EPRLAB - AMPROJ ) * ( EPRLAB + AMPROJ ) )
      PXTARG = PXCHCK - PPRLAB * TXX
      PYTARG = PYCHCK - PPRLAB * TYY
      PZTARG = PZCHCK - PPRLAB * TZZ
      GAM    = ETRLAB / AMTAR
      BGX    = PXTARG / AMTAR
      BGY    = PYTARG / AMTAR
      BGZ    = PZTARG / AMTAR
      PPHELP = ( BGX * TXX + BGY * TYY + BGZ * TZZ ) * PPRLAB
      ETAPCM = EPRLAB - PPHELP / ( GAM + ONEONE )
      PXPROJ = PPRLAB * TXX - BGX * ETAPCM
      PYPROJ = PPRLAB * TYY - BGY * ETAPCM
      PZPROJ = PPRLAB * TZZ - BGZ * ETAPCM
      UUOLD  = PXPROJ / PPROJX
      VVOLD  = PYPROJ / PPROJX
      WWOLD  = PZPROJ / PPROJX
      SINT02 = UUOLD**2 + VVOLD**2
      IF ( SINT02 .LE. ANGLSQ ) THEN
         LSMPAN = .TRUE.
         SINTH0 = ZERZER
         COSPH0 = ONEONE
         SINPH0 = ZERZER
      ELSE
         LSMPAN = .FALSE.
         SINTH0 = SQRT (SINT02)
         COSPH0 = UUOLD / SINTH0
         SINPH0 = VVOLD / SINTH0
      END IF
      PLABS  = PPROJX
      ELABS  = EPROJX
      IF ( LSEACL .OR. .NOT. LDIFFR (KPTOIP(KP)) .OR. PLABS .LE. PTHDFF
     &    ) THEN
         FRUUN(1) = ONEONE
      ELSE
         IF ( RN1GSC .GE. ZERZER ) THEN
            CALL GRNDM(FRUUN,1)
            IF ( FRUUN(1) .LT. HLFHLF ) THEN
               FRUUN(1) = RN1GSC
            ELSE
               FRUUN(1) = RN2GSC
            END IF
         ELSE
            CALL GRNDM(FRUUN,1)
         END IF
      END IF
 
      IF ( UMO * UMO .LT. URMIN2 ) THEN
         IF ( URMIN2 - UMO2 .LT. 0.1D+00 * URMIN2 .AND. IIBAR (KP) .LT.
     &      0 ) GO TO 1550
         IF ( .NOT. LSEACL .AND. PLABS .GT. 2.D+00 ) GO TO 1550
         WRITE ( LUNOUT,* )' Ferevv: Umo2 < Urmin2 !!',UMO*UMO,URMIN2,
     &                       KP,KT,AMPROJ,AMTAR
         WRITE ( LUNERR,* )' Ferevv: Umo2 < Urmin2 !!',UMO*UMO,URMIN2,
     &                       KP,KT,AMPROJ,AMTAR
         NHAD = 2
         NREH(1) = KT
         PXH(1) = 0.D+00
         PYH(1) = 0.D+00
         PZH(1) = 0.D+00
         HEPH(1) = AMTAR
         AMH (1) = AMTAR
         IBARH (1) = IIBAR (KT)
         ICHH  (1) = IICH  (KT)
         ANH   (1) = ANAME (KT)
         IF ( LSEACL .OR. IVFLAG .EQ. 0 ) THEN
            NREH(2) = 23
         ELSE
            NREH(2) = KP
         END IF
         PXH(2) = 0.D+00
         PYH(2) = 0.D+00
         PZH(2) = PLABS
         HEPH(2) = ELABS
         AMH (2) = AMPROJ
         IBARH (2) = IIBAR (NREH(2))
         ICHH  (2) = IICH  (NREH(2))
         ANH   (2) = ANAME (NREH(2))
         GO TO 1122
      END IF
 1550 CONTINUE
      IF ( FRUUN(1) .LE. FRDIFF ) THEN
         CALL DIFEVV ( NHAD, KP, KT, PLABS, ELABS, UMO )
         LEVDIF = .TRUE.
      ELSE IF ( .NOT. LSEACL .AND. PLABS .LT. DBLE(IJNUCR(KPTOIP(KP)))
     &        * 0.8D+00 * PTHRSH (KPTOIP(KP)) ) THEN
         IF ( KP .EQ. 26 ) THEN
            KPP = 23
         ELSE
            KPP = KP
         END IF
         CALL HEVHIN ( NHAD, KPP, KT, PLABS, ELABS, UMO )
         LHADRI = .TRUE.
      ELSE
         IF ( .NOT. LSEACL ) THEN
            LEVDIF = .FALSE.
            LHADRI = .FALSE.
         END IF
         CALL HADEVV ( NHAD, KP, KT, PLABS, ELABS, UMO )
      END IF
1122  CONTINUE
      DO 2000 I=1,NHAD
         IF ( LSMPAN ) THEN
            PZH (I) = WWOLD * PZH (I)
         ELSE
            PLRX = PXH (I) * COSPH0 * WWOLD - PYH (I) * SINPH0
     &           + PZH (I) * UUOLD
            PLRY = PXH (I) * SINPH0 * WWOLD + PYH (I) * COSPH0
     &           + PZH (I) * VVOLD
            PLRZ = - PXH (I) * SINTH0 + PZH (I) * WWOLD
            PXH (I) = PLRX
            PYH (I) = PLRY
            PZH (I) = PLRZ
         END IF
         CALL ALTRA ( GAM, BGX, BGY, BGZ, PXH(I), PYH(I), PZH(I),
     &                HEPH(I), PLR, PLRX, PLRY, PLRZ, ELR )
         PXH(I) = PLRX
         PYH(I) = PLRY
         PZH(I) = PLRZ
         HEPH(I) = ELR
2000  CONTINUE
 
      B1BAMJ = B1SAVE
      B2BAMJ = B2SAVE
      RETURN
      END
