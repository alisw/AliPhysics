*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*$ CREATE NUCEVV.FOR
*COPY NUCEVV
*
*=== nucevv ===========================================================*
*
      SUBROUTINE NUCEVV ( NNHAD, KPROJ, PPPROJ, EKPROJ, TXI, TYI, TZI )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*----------------------------------------------------------------------*
*     Nucevt90: new version by A. Ferrari INFN-Milan and CERN-SPS
*
*       Besides code cleaning, some changes have been made to extract
*       the quantities needed for a correct energy, momentum, electric
*       charge and baryonic charge conservation. This quantities are put
*       in /balanc/ common
*----------------------------------------------------------------------*
*
C
C**************************************************************
C        /NUCPAR/
C        PARTICLES CREATED IN NUCEVT
C           PXNU,PYNU,PZNU = X-, Y- AND Z-COMPONENTS OF THE
C           LAB MOMENTUM OF THE SECONDARY. COORDINATE SYSTEM
C           DEFINED BY THE PRIMARY PARTICLE.
C           HEPNU    =  TOTAL ENERGY IN THE LAB
C           AMNU     =  REST MASS
C           ICHNU    =  CHARGE
C           IBARNU   =  BARYONIC NUMBER
C           ANNU     =  NAME OF THE PARTICLE
C           NRENU    =  TYPE NUMBER OF THE PARTICLE
C**************************************************************
C
#include "geant321/balanc.inc"
#include "geant321/corinc.inc"
#include "geant321/depnuc.inc"
#include "geant321/hadpar.inc"
#include "geant321/inpdat2.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucpar.inc"
#include "geant321/parevt.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
      COMMON /FKPRIN/ IPRI, INIT
*  The following dimension statement is now obsolete
*     DIMENSION XQT(20), XDQT(20), IFQT(3,20)
      REAL RNDM(1)
      LOGICAL LDSAVE
      DIMENSION PTHRSH (39)
      SAVE PTHRSH
      DATA PTHRSH / 16*5.D+00,2*2.5D+00,5.D+00,3*2.5D+00,8*5.D+00,
     &              9*2.5D+00 /
*
*
*
*----------------------------------------------------------------------*
*         Eproj  = total energy of the projectile
*         Amproj = mass energy of the projectile
*         Ekproj = kinetic energy of the projectile
*         Pproj  = momentum of the projectile
*         Eeproj = original total energy of the projectile
*         Ppproj = original momentum of the projectile
*         Ibproj = barionic charge of the projectile
*         Icproj = electric charge of the projectile
*         Nnhad  = number of particles produced in Nucevt
*         Nevt   = number of high energy collisions
*         Atemp  = local variable with the mass number of the residual
*                  nucleus : it is updated at any interaction and alrea-
*                  dy includes contributions from cascade particles
*                  (even though actually they are produced after the
*                   high energy interactions: since there is no real
*                   correlation except for the one with Nsea, this does
*                   not represent a problem)
*         Ztemp  = same as Atemp for the atomic number
*----------------------------------------------------------------------*
*
      TXX = TXI
      TYY = TYI
      TZZ = TZI
      PPROJ = PPPROJ
      AMPROJ= AM (KPROJ)
      EPROJ = EKPROJ + AMPROJ
      EEPROJ= EPROJ
C
      IBPROJ= IBAR (KPROJ)
      ICPROJ= ICH  (KPROJ)
*  Set Atemp and Ztemp to their initial values
      ATEMP = IBTAR  - IGREYP - IGREYN
      ZTEMP = ICHTAR - IGREYP
      NNHAD = 0
      NEVT  = 0
* Now Nsea is sampled inside Corrin and passed through common /corrinc
      IF ( NSEA .EQ. 0 ) GO TO 1000
      IF ( INIT .EQ. 1 ) WRITE(LUNOUT,1)NNHAD,KPROJ,PPROJ,
     *                   AMPROJ,EPROJ,ANUAV,NSEA
   1  FORMAT (1X,2I5,4F12.5,I10)
*
*  We have now sampled Nsea, the number of quark-antiquark pairs
*  interacting besides the valence interaction (total interactions
*  = Nsea+1
*
*or      IF (INIT.EQ.1) WRITE(LUNOUT,4)XO,(XSEA(I),XASEA(I),I=1,NSEA)
*or   4  FORMAT (10F12.6)
*  +-------------------------------------------------------------------*
*  |  Sample effective meson nucleus collisions
*  |  Em, im, ekm, pm refer to the meson quantities
      DO 6 I = 1, NSEA
         EM = EKPROJ * ( XSEA(I) + XASEA(I) )
*  |  +----------------------------------------------------------------*
*  |  |  Selection of the ddbar or qqbar composition
         CALL GRNDM(RNDM,1)
         IF ( RNDM (1) .LT. 0.5D+00 ) THEN
            IM = 23
*  |  |
*  |  +----------------------------------------------------------------*
*  |  |
         ELSE
            IM = 26
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         AMPIO  = AM (IM)
         EKM = EM - AMPIO
*  |  +----------------------------------------------------------------*
*  |  |  Energy and/or momentum is too low!!!
         IF ( EKM .LT. ETHSEA .OR. PPROJ .LT. PTHRSH (KPTOIP(KPROJ)) )
     &      THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( I .LT. NSEA ) THEN
               II = I + 1
               XSEA(II) = XSEA(II)  + XSEA(I)
               XASEA(II)= XASEA(II) + XASEA(I)
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            KT = IJTARG (I)
*  |  |  +-------------------------------------------------------------*
*  |  |  | Kt is the index of the target nucleon (1=proton,8=neutron)
            IF ( KT .EQ. 1 ) THEN
               ZNOW  = ZNOW + 1.D+00
               KTARP = KTARP - 1
               ZNCOLL = ZNCOLL - 1.D+00
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               KTARN = KTARN - 1
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            ANOW   = ANOW   + 1.D+00
            ANCOLL = ANCOLL - 1.D+00
            GO TO 6
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         AMPIO2 = AMPIO * AMPIO
         EPROLD = EKPROJ + AMPROJ
         PPROLD = PPROJ
*  | After selection of the meson energy, Ekproj is updated
         EKPROJ = EKPROJ - EM
         EPROJ  = EKPROJ + AMPROJ
*  |  +----------------------------------------------------------------*
*  |  | All ekproj energy is transferred to em if ekproj < 4 GeV
*  |  | Now modified since it causes troubles to energy conservation by
*  |  | A. Ferrari: force them to have only a valence interaction
         IF ( EKPROJ .LT. 3.99D0 ) THEN
            EKPROJ = EKPROJ + EM
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            DO 666 IN = I, NSEA
               KT = IJTARG (IN)
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  | Kt is the index of the target nucleon (1=proton,8=neutron)
               IF ( KT .EQ. 1 ) THEN
                  ZNOW  = ZNOW + 1.D0
                  KTARP = KTARP - 1
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               ELSE
                  KTARN = KTARN - 1
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               ANOW  = ANOW + 1.D0
 666        CONTINUE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            GO TO 66
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Now the wounded nucleus is selected inside Corevt
         KT    = IJTARG (I)
         ATEMP = ATEMP - 1.D+00
         IF ( KT .EQ. 1 ) ZTEMP = ZTEMP - 1.D+00
*  |  +---------------------------------------------------------------*
*  |  |
         IF ( NINT ( ATEMP ) .EQ. 1 ) THEN
            LLASTN = .TRUE.
            KTLAST = KT
            KTINC  = IJTARG ( NSEA + 1 )
            PM = SQRT ( EM**2 - AMPIO2 )
            AMLAST = AM (KTLAST)
            AMINC  = AM (KTINC)
            UMIN2  = ( AMINC + AMLAST )**2
            ITJ = MIN ( KTINC, 2 )
            DELTAE = V0WELL (ITJ) - EFRMAV (ITJ)
*  |  |  +------------------------------------------------------------*
*  |  |  |   Selection of the invariant mass of the two last nucleon
*  |  |  |   system
2020        CONTINUE
               EKPROJ = EKPROJ - DELTAE
               EFRM   = EFRM  - DELTAE
               ELEFT  = ETTOT - EINTR - EUZ - EKPROJ - AMPROJ - EM
               PPROJ  = SQRT ( EKPROJ * ( EKPROJ + 2.D+00 * AMPROJ ) )
               PXLEFT = PXTTOT - PXINTR - PUX - ( PPROJ + PM ) * TXX
               PYLEFT = PYTTOT - PYINTR - PUY - ( PPROJ + PM ) * TYY
               PZLEFT = PZTTOT - PZINTR - PUZ - ( PPROJ + PM ) * TZZ
               UMO2 = ELEFT**2 - PXLEFT**2 - PYLEFT**2 - PZLEFT**2
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( UMO2 .LE. UMIN2 ) THEN
                  PPDTPL = PXLEFT * TXX + PYLEFT * TYY + PZLEFT * TZZ
                  DELTAE = 0.51D+00 * ( UMIN2 - UMO2 ) / ( ELEFT -
     &                     PPDTPL * ( EKPROJ + AMPROJ ) / PPROJ )
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |   Skip the sea interaction if deltae < 0
                  IF ( DELTAE .LT. 0.D+00 ) THEN
                     EKPROJ = EKPROJ + EM
                     ATEMP = ATEMP + 1.D+00
                     IF ( KT .EQ. 1 ) ZTEMP = ZTEMP + 1.D+00
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  |
                     DO 777 IN = I, NSEA
                        KT = IJTARG (IN)
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |  Kt is the index of the target nucleon
*  |  |  |  |  |  |  | (1=proton,8=neutron)
                        IF ( KT .EQ. 1 ) THEN
                           ZNOW  = ZNOW + 1.D0
                           KTARP = KTARP - 1
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
*  |  |  |  |  |  |  |
                        ELSE
                           KTARN = KTARN - 1
                        END IF
*  |  |  |  |  |  |  |
*  |  |  |  |  |  |  +-------------------------------------------------*
                        ANOW  = ANOW + 1.D0
 777                 CONTINUE
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
                     PPROJ  = SQRT ( EKPROJ * ( EKPROJ + 2.D+00
     &                      * AMPROJ ) )
                     GO TO 66
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  ELSE IF ( DELTAE .GE. EKPROJ ) THEN
                     WRITE ( LUNERR,* )' Nucevv: sea call impossible ',
     &                                 ' to get Umin2, go to resampling'
                     LRESMP = .TRUE.
                     RETURN
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
                  GO TO 2020
*  |  |  |
*  |  |  +-<|--<--<--<--< We need more invariant mass!!
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Now we will divide Eleft and Pleft between the two
*  |  |  left nucleons!
            UMO = SQRT ( UMO2 )
            ELACMS = 0.5D+00 * ( UMO2 + AMLAST**2 - AMINC**2 )
     &             / UMO
            EINCMS = UMO - ELACMS
            PCMS   = SQRT ( ELACMS**2 - AMLAST**2 )
            GAMCM = ELEFT  / UMO
            ETAX  = PXLEFT / UMO
            ETAY  = PYLEFT / UMO
            ETAZ  = PZLEFT / UMO
            CALL RACO ( CXXINC, CYYINC, CZZINC )
            PCMSX = PCMS * CXXINC
            PCMSY = PCMS * CYYINC
            PCMSZ = PCMS * CZZINC
*  |  |  Now go back from the CMS frame to the lab frame!!!
*  |  |  First the "inc" nucleon:
            ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
            EKINC  = GAMCM * EINCMS + ETAPCM - AMINC
            PHELP  = ETAPCM / (GAMCM + 1.D+00) + EINCMS
            PXXINC = PCMSX + ETAX * PHELP
            PYYINC = PCMSY + ETAY * PHELP
            PZZINC = PCMSZ + ETAZ * PHELP
*  |  |  Now the "last" nucleon
            EKLAST = GAMCM * ELACMS - ETAPCM - AMLAST
            PHELP  = - ETAPCM / (GAMCM + 1.D+00) + ELACMS
            PXLAST = - PCMSX + ETAX * PHELP
            PYLAST = - PCMSY + ETAY * PHELP
            PZLAST = - PCMSZ + ETAZ * PHELP
            TVEUZ  = 0.D+00
            TVGRE0 = 0.D+00
            EINCT  = EINCP + EINCN
            IF ( EINCT .GT. 0.D+00 ) THEN
               EINCP  = EINCP - EINCP / EINCT * TVGREY
               EINCN  = EINCN - EINCN / EINCT * TVGREY
            END IF
            TVGREY = 0.D+00
            PSEA = PM + PPROJ - PPROLD
            LDSAVE = LDIFFR (KPTOIP(IM))
            LDIFFR (KPTOIP(IM)) = .FALSE.
*  |  |  The last argument for ferevv is set to 0 to signal that it
*  |  |  is not a true valence call
            IVFLAG = 0
            CALL FEREVV ( NHAD, IM, KT, PM, EKM, TXX, TYY, TZZ, ATEMP,
     &                    ZTEMP, IVFLAG )
            LDIFFR (KPTOIP(IM)) = LDSAVE
            IF ( LRESMP ) RETURN
*  |  |
*  |  +---------------------------------------------------------------*
*  |  |
         ELSE
            TXXOLD = TXX
            TYYOLD = TYY
            TZZOLD = TZZ
            PPROJ  = SQRT ( EKPROJ * ( EKPROJ + 2.D+00 * AMPROJ ) )
            IF ( PPROJ .GT. PPROLD ) THEN
               PPROJ = PPROLD * EPROJ / EPROLD
            END IF
            CALL FERSEA ( NHAD, IM, KT, EKM, TXX, TYY, TZZ, ATEMP,
     &                    ZTEMP, KPROJ, EPROJ, PPROJ, EPROLD, PPROLD )
            IF ( LRESMP ) RETURN
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         NEVT = NEVT + 1
*  |  +----------------------------------------------------------------*
*  |  |
         IF ( NNHAD + NHAD .GT. MXPNUC ) THEN
            WRITE (LUNOUT,1101) NNHAD
1101        FORMAT('  NHAD IN NUCEVT TOO BIG CHANGE DIMENSION',
     *             '  NNHAD=',I10)
            STOP
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
 
*  |  +----------------------------------------------------------------*
*  |  |          Looping over the produced particles!!!
         DO 7 JJ = 1, NHAD
            NNHAD = NNHAD + 1
            J     = NNHAD
            PXNU(J)  = PXH(JJ)
            PYNU(J)  = PYH(JJ)
            PZNU(J)  = PZH(JJ)
            HEPNU(J) = HEPH(JJ)
            AMNU(J)  = AMH(JJ)
            IBARNU(J)= IBARH(JJ)
            ICHNU(J) = ICHH(JJ)
            NRENU(J) = NREH(JJ)
            ANNU(J)  = ANH(JJ)
*  |  |          Updating energy and momentum accumulators
            PUX = PUX+PXNU(J)
            PUY = PUY+PYNU(J)
            PUZ = PUZ+PZNU(J)
            EUZ = EUZ+HEPNU(J)
            ICU = ICU+ICHNU(J)
            IBU = IBU+IBARNU(J)
   7     CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
   6  CONTINUE
*  |
*  +-------------------------------------------------------------------*
  66  CONTINUE
*
*  Computing the actual pproj
*
      EPROJ = EKPROJ + AMPROJ
1000  CONTINUE
*  +-------------------------------------------------------------------+
*  |  Now modified because of troubles to energy conservation by
*  |  A. Ferrari: if the incident projectile is a proton or a neutron
*  |  maybe it was stopped in the nucleus, else it is added to the stack,
*  |  other particles are also added to the secondary stack.
*  |  May be is not physical, but better than nothing!!! A different
*  |  solution maybe to let meson and other baryons decay at rest or
*  |  also force them to have only a valence interaction or also
*  |  convert a charged meson plus a residual neutron/proton to a
*  |  proton/neutron and add the extra energy to the sea meson and
*  |  so on        ????????????????????
*  |  Now the condition Pproj .ge. .4 is always fulfilled!!!!!
      IF (PPROJ .LT. 4.D0) THEN
*  |  +---------------------------------------------------------------*
*  |  |
         DKPR30 = KPROJ-30
         IF ( PPROJ - 1.5D+00 * SIGN ( ONEONE, DKPR30 )
     &        .LT. 4.D0 ) THEN
            WRITE ( LUNERR,* )' Nucevt: Pproj < 4 GeV/c!!!', PPROJ,
     &                          KPROJ
         END IF
*  |  |
*  |  +---------------------------------------------------------------*
      END IF
*  |
*  +------------------------------------------------------------------*
*   Now the wounded nucleus is selected inside Corevt
      KT = IJTARG ( NSEA + 1 )
*   Update the Nsea value
      NSEA = NEVT
*or      IF (INIT.EQ.1) WRITE(LUNOUT,162)NEVT,NHAD,IM,KPROJ,KT,PM,
*or     &EM,EKM,PM,EKPROJ,PPROJ,EPROJ
*
*  Now ferevt performs the interaction: no longer a meson but the actual
*  projectile
*
      ATEMP = ATEMP - 1.D+00
      IF ( KT .EQ. 1 ) ZTEMP = ZTEMP - 1.D+00
*  +-------------------------------------------------------------------*
*  |
      IF ( NINT ( ATEMP ) .EQ. 1 ) THEN
         LLASTN = .TRUE.
         KTLAST = KT
         IF ( ZNOW .GT. 0.D+00 ) THEN
            KTINC  = 1
         ELSE
            KTINC  = 8
         END IF
         AMLAST = AM (KTLAST)
         AMINC  = AM (KTINC)
         UMIN2  = ( AMINC + AMLAST )**2
         ITJ = MIN ( KTINC, 2 )
*  |  +----------------------------------------------------------------*
*  |  |   Selection of the invariant mass of the two last nucleon
*  |  |   system
2040     CONTINUE
            ELEFT  = ETTOT - EINTR - EUZ - EKPROJ - AMPROJ
            PXLEFT = PXTTOT - PXINTR - PUX - PPROJ * TXX
            PYLEFT = PYTTOT - PYINTR - PUY - PPROJ * TYY
            PZLEFT = PZTTOT - PZINTR - PUZ - PPROJ * TZZ
            UMO2 = ELEFT**2 - PXLEFT**2 - PYLEFT**2 - PZLEFT**2
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( UMO2 .LE. UMIN2 ) THEN
               PPDTPL = PXLEFT * TXX + PYLEFT * TYY + PZLEFT * TZZ
               DELTAE = 0.51D+00 * ( UMIN2 - UMO2 ) / ( ELEFT -
     &                  PPDTPL * ( EKPROJ + AMPROJ ) / PPROJ )
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               IF ( DELTAE .GE. EKPROJ .OR. DELTAE .LT. 0.D+00 ) THEN
                  WRITE ( LUNERR,* )' Nucevv: valence call impossible ',
     &                              ' to get Umin2, go to resampling'
                  LRESMP = .TRUE.
                  RETURN
               END IF
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               EKPROJ = EKPROJ - DELTAE
               PLA    = PPROJ
               PPROJ  = SQRT ( EKPROJ * ( EKPROJ + 2.D+00 * AMPROJ ) )
               PLA    = PPROJ - PLA
               EFRM   = EFRM  - DELTAE
               PXFRM  = PXFRM + TXX * PLA
               PYFRM  = PYFRM + TYY * PLA
               PZFRM  = PZFRM + TZZ * PLA
               GO TO 2040
*  |  |
*  |  +-<|--<--<--<--< We need more invariant mass!!
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |
*  |  +----------------------------------------------------------------*
*  |  Now we will divide Eleft and Pleft between the two
*  |  left nucleons!
         UMO = SQRT ( UMO2 )
         ELACMS = 0.5D+00 * ( UMO2 + AMLAST**2 - AMINC**2 ) / UMO
         EINCMS = UMO - ELACMS
         PCMS   = SQRT ( ELACMS**2 - AMLAST**2 )
         GAMCM = ELEFT  / UMO
         ETAX  = PXLEFT / UMO
         ETAY  = PYLEFT / UMO
         ETAZ  = PZLEFT / UMO
         CALL RACO ( CXXINC, CYYINC, CZZINC )
         PCMSX = PCMS * CXXINC
         PCMSY = PCMS * CYYINC
         PCMSZ = PCMS * CZZINC
*  |  Now go back from the CMS frame to the lab frame!!!
*  |  First the "inc" nucleon:
         ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
         EKINC  = GAMCM * EINCMS + ETAPCM - AMINC
         PHELP  = ETAPCM / (GAMCM + 1.D+00) + EINCMS
         PXXINC = PCMSX + ETAX * PHELP
         PYYINC = PCMSY + ETAY * PHELP
         PZZINC = PCMSZ + ETAZ * PHELP
*  |  Now the "last" nucleon
         EKLAST = GAMCM * ELACMS - ETAPCM - AMLAST
         PHELP  = - ETAPCM / (GAMCM + 1.D+00) + ELACMS
         PXLAST = - PCMSX + ETAX * PHELP
         PYLAST = - PCMSY + ETAY * PHELP
         PZLAST = - PCMSZ + ETAZ * PHELP
         TVEUZ  = 0.D+00
         TVGRE0 = 0.D+00
         EINCT  = EINCP + EINCN
         IF ( EINCT .GT. 0.D+00 ) THEN
            EINCP  = EINCP - EINCP / EINCT * TVGREY
            EINCN  = EINCN - EINCN / EINCT * TVGREY
         END IF
         TVGREY = 0.D+00
      END IF
*  |
*  +-------------------------------------------------------------------*
*  The last argument for ferevv is set to 1 to signal that it
*  is a true valence call
      IVFLAG = 1
      CALL FEREVV ( NHAD, KPROJ, KT, PPROJ, EKPROJ, TXX, TYY, TZZ,
     &              ATEMP, ZTEMP, IVFLAG )
      IF ( LRESMP ) RETURN
      NEVT = NEVT + 1
*  +-------------------------------------------------------------------*
*  |
      IF ( NNHAD + NHAD .GT. MXPNUC ) THEN
         WRITE (LUNOUT,1101) NNHAD
         STOP
      END IF
*  |
*  +-------------------------------------------------------------------*
 
*  +-------------------------------------------------------------------*
*  |       Looping over the produced particles
      DO 8 JJ=1,NHAD
         NNHAD = NNHAD+1
         J = NNHAD
         PXNU(J)  = PXH(JJ)
         PYNU(J)  = PYH(JJ)
         PZNU(J)  = PZH(JJ)
         HEPNU(J) = HEPH(JJ)
         AMNU(J)  = AMH(JJ)
         IBARNU(J)= IBARH(JJ)
         ICHNU(J) = ICHH(JJ)
         NRENU(J) = NREH(JJ)
         ANNU(J)  = ANH(JJ)
*  |             Updating energy and momentum accumulators
         PUX = PUX+PXNU(J)
         PUY = PUY+PYNU(J)
         PUZ = PUZ+PZNU(J)
         EUZ = EUZ+HEPNU(J)
         ICU = ICU+ICHNU(J)
         IBU = IBU+IBARNU(J)
   8  CONTINUE
*  |
*  +-------------------------------------------------------------------*
 888  CONTINUE
*     IF (PPROJ .LT. 4.D0) INIT=1
**** print and test energy conservation
*
      ETOT  = EEPROJ + KTARP * ( AMNUCL (1) - EBNDNG (1) )
     &      + KTARN * ( AMNUCL (2) - EBNDNG (2) ) + EFRM
      ICTOT = ICH  (KPROJ) + KTARP
      IBTOT = IBAR (KPROJ) + KTARP + KTARN
      EMIN  = 1.D-10 * ETOT
      PMIN  = 1.D-10 * PPPROJ
*  +-------------------------------------------------------------------*
*  |
      IF (ICTOT .NE. ICU .OR. IBTOT .NE. IBU) THEN
         LRESMP = .TRUE.
         WRITE(LUNERR,*)' NUCEVT-FEREVT CHARGE CONSERVATION FAILURE:',
     &                  ' ICU,ICTOT,IBU,IBTOT',ICU,ICTOT,IBU,IBTOT
      END IF
*  |
*  +-------------------------------------------------------------------*
      PXCHCK = PPPROJ * TXI + PSEA * TXX
      PYCHCK = PPPROJ * TYI + PSEA * TYY
      PZCHCK = PPPROJ * TZI + PSEA * TZZ
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS ( PUX - PXFRM - PXCHCK ) .GT. PMIN ) THEN
         LRESMP = .TRUE.
         WRITE(LUNERR,*)' NUCEVT-FEREVT PX CONSERVATION FAILURE:',
     &                  ' PUX,PXFRM+PXCHCK',PUX,PXFRM+PXCHCK
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS ( PUY - PYFRM - PYCHCK ) .GT. PMIN ) THEN
         WRITE(LUNERR,*)' NUCEVT-FEREVT PY CONSERVATION FAILURE:',
     &                  ' PUY,PYFRM+PYCHCK*TYY',PUY,PYFRM+PYCHCK
         LRESMP = .TRUE.
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS ( PUZ - PZFRM - PZCHCK ) .GT. PMIN ) THEN
         LRESMP = .TRUE.
         WRITE(LUNERR,*)' NUCEVT-FEREVT PZ CONSERVATION FAILURE:',
     &                  ' PUZ,PZFRM+PZCHCK',PUZ,PZFRM+PZCHCK
      END IF
*  |
*  +-------------------------------------------------------------------*
*  +-------------------------------------------------------------------*
*  |
      IF ( ABS ( ETOT - EUZ ) .GT. EMIN ) THEN
         LRESMP = .TRUE.
         WRITE(LUNERR,*)' NUCEVT-FEREVT E CONSERVATION FAILURE:',
     &                  ' EUZ,ETOT',EUZ,ETOT
      END IF
*  |
*  +-------------------------------------------------------------------*
      EUZ0 = EUZ-NEVT*AM(1)
      IF (INIT.NE.1) GO TO 11
      WRITE(LUNOUT,12)NNHAD,KPROJ,PPPROJ,EPROJ,PUX,PUY,PUZ,EUZ0,
     *ICU,IBU,NNHAD
  12  FORMAT (1X,2I5,6F12.6,3I5)
*  +-------------------------------------------------------------------*
*  |
      DO 13 I=1,NNHAD
         WRITE(LUNOUT,14)I,NRENU(I),ICHNU(I),IBARNU(I),ANNU(I),PXNU(I),
     *                   PYNU(I),PZNU(I),HEPNU(I),AMNU(I)
  14     FORMAT (1X,4I5,A8,6F12.6)
  13  CONTINUE
*  |
*  +-------------------------------------------------------------------*
      INIT=0
  11  CONTINUE
      RETURN
      END
