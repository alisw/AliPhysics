*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:05  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*=== dres =============================================================*
*                                                                      *
      SUBROUTINE FKDRES ( M2, M3, T1, U, EREC, LOPPAR, JFISS )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*  New version of DRES created  by A.Ferrari & P.Sala, INFN - Milan    *
*                                                                      *
*  Last change  on  10-apr-93  by  Alfredo Ferrari, INFN - Milan       *
*                                                                      *
*  Dres93: Dres91 plus the RAL fission model taken from LAHET thanks   *
*          to R.E.Prael                                                *
*  Dres91: new version from A. Ferrari and P. Sala, INFN - Milan       *
*          This routine has been adapted from the original one of the  *
*          Evap-5 module (KFA - Julich). Main modifications concern    *
*          with kinematics which is now fully relativistic and with    *
*          the treatment of few nucleons nuclei, which are now frag-   *
*          mented, even though in a very rough manner. Changes have    *
*          been made also to other routines of the Evap-5 package      *
*                                                                      *
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
*                                                                      *
*  Input variables:                                                    *
*     M2 = Mass number of the residual nucleus                         *
*     M3 = Atomic number of the residual nucleus                       *
*     T1 = Excitation energy of the residual nucleus before evaporation*
*     U  = Excitation energy of the residual nucleus after evaporation *
*     Erec = Recoil kinetic energy of the residual nucleus             *
*            The recoil direction is given by Coslbr (i)               *
*                                                                      *
*  Significant variables:                                              *
*     JA = Present mass number of the residual nucleus                 *
*     JZ = Present atomic number of the residual nucleus               *
*     Smom1 = Energy accumulators for the six types of evaporated      *
*             particles                                                *
*                                                                      *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*    !!!! Please note that the following variables concerning !!!!     *
*    !!!! with the present residual nucleus must be set before!!!!     *
*    !!!! entering DRES91: Ammres, Ptres                      !!!!     *
*    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     *
*                                                                      *
*----------------------------------------------------------------------*
*
C---------------------------------------------------------------------
C SUBNAME = DRES --- EVAPORATION
C           EVAPORATION DATA SHOUD BE READ ON INPUT STAGE
C---------------------------------------------------------------------
C*****A LA EVAP III(TWA,8-68)
#include "geant321/eva0.inc"
#include "geant321/eva1.inc"
#include "geant321/forcn.inc"
#include "geant321/higfis.inc"
#include "geant321/inpflg.inc"
#include "geant321/labcos.inc"
#include "geant321/nucdat.inc"
#include "geant321/resnuc.inc"
      DIMENSION ZMASS (6), Z2MASS(6), C(3), Q(0:6), FLKCOU(6), CCOUL(6),
     &          THRESH(6), SMALLA(6), R(6), S  (6), SOS   (6), STRUN(6),
     &          EYE1  (6), EYE0  (6), SMOM1    (6), BNMASS(6)
      DIMENSION CORRRR(6)
      REAL RNDM(2)
      LOGICAL LOPPAR, PENBAR, LFIRST
      PARAMETER (EXPMIN=-88,EXPMAX=88)
      SAVE ZMASS, Z2MASS, EMHN, EMNUM, UM, AMUMEV, AMEMEV, QBRBE8,
     &     BNMASS, IEVEVP, NBE8, NRNEEP, LFIRST
      DATA IEVEVP / 0 /
      DATA LFIRST / .TRUE. /
*
      IEVEVP = IEVEVP + 1
C-------------------------------------- 1.ST CALL INIT
      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         FKEY   = ZERZER
         NBE8   = 0
         NRNEEP = 0
         EXMASS(1) = 1.D+03 * ( AMNEUT - AMUAMU )
         EXMASS(2) = FKENER ( ONEONE, ONEONE )
         EXMASS(3) = FKENER ( TWOTWO, ONEONE )
         EXMASS(4) = FKENER ( THRTHR, ONEONE )
         EXMASS(5) = FKENER ( THRTHR, TWOTWO )
         EXMASS(6) = FKENER ( FOUFOU, TWOTWO )
         ZMASS(1) = 1.D+03 * AMUAMU + EXMASS (1)
         ZMASS(2) = 1.D+03 * AMUAMU + EXMASS (2)
         ZMASS(3) = 2.D+03 * AMUAMU + EXMASS (3)
         ZMASS(4) = 3.D+03 * AMUAMU + EXMASS (4)
         ZMASS(5) = 3.D+03 * AMUAMU + EXMASS (5)
         ZMASS(6) = 4.D+03 * AMUAMU + EXMASS (6)
         BNMASS (1) = 0.D+00
         BNMASS (2) = 0.D+00
         BNMASS (3) = ZMASS (1) + ZMASS (2) - ZMASS (3)
         BNMASS (4) = TWOTWO * ZMASS (1) + ZMASS (2) - ZMASS (4)
         BNMASS (5) = ZMASS (1) + TWOTWO * ZMASS (2) - ZMASS (5)
         BNMASS (6) = TWOTWO * ( ZMASS (1) + ZMASS (2) ) - ZMASS (6)
         DO 1234 KK = 1,6
            Z2MASS (KK) = ZMASS (KK) * ZMASS (KK)
1234     CONTINUE
         AMUMEV = 1.D+03 * AMUAMU
         AMEMEV = 1.D+03 * AMELEC
         QBRBE8 = FKENER ( EIGEIG, FOUFOU ) - TWOTWO * EXMASS (6)
         EMN = 1.D+03 * AMNEUT
         EMH = ZMASS (2)
         TMP16 = 16.D+00
         UM  = AMUMEV + FKENER ( TMP16, EIGEIG ) / 16.D+00
         EMHN  = EMH - EMN
         EMNUM = EMN - UM
      END IF
*  |  End of initialization:
*  +-------------------------------------------------------------------*
C     --------------------------------- START OF PROCESS
*  +-------------------------------------------------------------------*
*  |  Initialize Npart and Smom if nothing has been already evaporated
*  |  for this event
      IF ( JFISS .LE. 0 ) THEN
         DO 775 I=1,6
            NPART(I) = 0
            SMOM1(I) = ZERZER
  775    CONTINUE
      END IF
*  |
*  +-------------------------------------------------------------------*
      JA = M2
      JZ = M3
      U  = T1
      RNMASS = 1.D+03 * AMMRES + U
* P2res and  Ptres are the squared momentum and the momentum of the
* residual nucleus (now in relativistic kinematics), Umo the
* invariant mass of the system!
      UMO  = RNMASS
      UMO2 = UMO * UMO
      ELBTOT  = RNMASS + EREC
      GAMCM   = ELBTOT / RNMASS
      ETACM   = 1.D+03 * PTRES / RNMASS
      HEVSUM  = ZERZER
 1000 CONTINUE
      LOPPAR = .FALSE.
*  +-------------------------------------------------------------------*
*  |                Check for starting data inconsistencies
      IF (JA-JZ .LT. 0) THEN
         WRITE(LUNOUT,6401)
         WRITE(LUNERR,6401)
 6401    FORMAT('1 Dres: cascade residual nucleus has mass no. less',
     &       ' than Z!!')
         RETURN
*  |
*  +-------------------------------------------------------------------*
*  |                Rough treatment for very few nucleon residual
*  |                nuclei. The basic ideas are:
*  |        a) as many as possible alpha particles are emitted
*  |        b) particles are emitted one per time leaving a residual
*  |           excitation energy proportional to number of nucleons
*  |           left in the residual nucleus (so we deal only with
*  |           two body kinematics)
*  |       T A K E   I N T O   A C C O U N T   T H A T   T H I S
*  |       T R E A T M E N T   I S   E X T R E M E L Y   R O U G H
*  |       T H E   T A S K   B E I N G   O N L Y   T O   S U P P L Y
*  |       S O M E T H I N G   T O   S H A R E   E N E R G Y   A N D
*  |       M O M E N T U M   A M O N G   A   F E W   F R A G M E N T S
      ELSE IF ( JA .LE. 6 .OR. JZ .LE. 2 ) THEN
*  | 1000 continue moved above according to FCA suggestion
*1000    CONTINUE
         JRESID = 0
         IF ( JA .GT. 4 ) GO TO 2000
*  |  +----------------------------------------------------------------*
*  |  | First check we are not concerning with a couple of neutrons or
*  |  | protons
         IF ( JA .EQ. 2 .AND. JZ .NE. 1 ) THEN
            JEMISS = 1 + JZ / 2
            JRESID = JEMISS
            RNMASS = ZMASS (JRESID)
            U = 0.D+00
            DELTU = UMO - 2.D+00 * ZMASS (JEMISS)
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( DELTU .LE. 0.D+00 ) THEN
               IF ( DELTU .LT. - 2.D+00 * ANGLGB * UMO ) THEN
                  WRITE ( LUNERR, * )' *** Dres: insufficient Umo for',
     &                               ' a nucleon couple', UMO,
     &                                 2.D+00 * ZMASS (JEMISS)
               END IF
               UMO = ( UMO + DELTU ) * ( 1.D+00 + ANGLGB )
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
            GO TO 2500
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  | Then check we are not concerning with one of the six
*  |  | standard particles
         DO 1700 J = 6, 1, -1
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( JZ .EQ. IZ (J) .AND. JA .EQ. IA (J) ) THEN
               HEVSUM = SMOM1(3) + SMOM1(5) + SMOM1(6) + SMOM1(4)
               GO TO ( 1100, 1100, 1600, 1500, 1400, 1300 ), J
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Proton or neutron, nothing can be done
 1100          CONTINUE
                  RETURN
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Alpha:
 1300          CONTINUE
                  DEUDEU = MAX ( ZERZER, U + TWOTWO * BNMASS (3)
     &                           - BNMASS (6) )
                  PROTRI = MAX ( ZERZER, U + BNMASS (4) - BNMASS (6) )
                  UEU3HE = MAX ( ZERZER, U + BNMASS (5) - BNMASS (6) )
                  QNORM  = DEUDEU + PROTRI + UEU3HE
*  |  |  |  |  If we cannot split then return
                  IF ( QNORM .LE. ZERZER ) RETURN
                  CALL GRNDM(RNDM,1)
                  V = RNDM (1)
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split or into two deuterons or a triton and a proton
*  |  |  |  |  | or a 3-He and a neutron: no account is made for
*  |  |  |  |  | Coulomb effects, probability is simply assumed
*  |  |  |  |  | proportional to reaction Qs
                  IF ( V .LT. DEUDEU / QNORM ) THEN
*  |  |  |  |  | Two deuterons selected
                     JEMISS = 3
                     JRESID = 3
                     RNMASS = ZMASS (3)
                     U = ZERZER
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split into a triton and a proton
                  ELSE IF ( V .LT. ( DEUDEU + PROTRI ) / QNORM ) THEN
                     JEMISS = 2
                     JRESID = 4
                     RNMASS = ZMASS (4)
                     U = ZERZER
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split into a 3-He and a neutron
                  ELSE
                     JEMISS = 1
                     JRESID = 5
                     RNMASS = ZMASS (5)
                     U = ZERZER
                     GO TO 2500
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  3-He:
 1400          CONTINUE
                  DEUPRO = MAX ( ZERZER, U + BNMASS (3) - BNMASS (5) )
                  PRPRNE = MAX ( ZERZER, U - BNMASS (5) )
                  QNORM  = DEUPRO + PRPRNE
*  |  |  |  |  If we cannot split then return
                  IF ( QNORM .LE. ZERZER ) RETURN
                  CALL GRNDM(RNDM,1)
                  V = RNDM (1)
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split or into a deuteron and a proton
*  |  |  |  |  | or into two protons and one neutron: no account is
*  |  |  |  |  | made for Coulomb effects, probability is simply assumed
*  |  |  |  |  | prportional to reaction Qs
                  IF ( V .LT. DEUPRO / QNORM ) THEN
*  |  |  |  |  | A deuteron and a proton selected
                     JEMISS = 2
                     JRESID = 3
                     RNMASS = ZMASS (3)
                     U = ZERZER
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split into 2 protons and 1 neutron: part of the exci-
*  |  |  |  |  | tation energy is conserved to allow the further
*  |  |  |  |  | splitting of the deuteron
                  ELSE
                     JEMISS = 2
                     JRESID = 0
                     FACT = ONEONE
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  | Loop to compute the residual excitation energy
 1450                CONTINUE
                        FACT = FACT * 0.6666666666666667D+00
*  |  |  |  |  |  | Erncm, Eepcm are the total energies of the residual
*  |  |  |  |  |  | nucleus and of the emitted particle in the CMS frame
                        U      = FACT * PRPRNE + BNMASS (3)
                        RNMASS = ZMASS (3) + U
                        ERNCM  = HLFHLF * ( UMO2 + RNMASS**2
     &                         - Z2MASS (JEMISS) ) / UMO
                        EEPCM  = UMO - ERNCM
                     IF ( EEPCM .LE. ZMASS (JEMISS) ) GO TO 1450
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
                     GO TO 2600
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Triton:
 1500          CONTINUE
                  DEUNEU = MAX ( ZERZER, U + BNMASS (3) - BNMASS (4) )
                  PRNENE = MAX ( ZERZER, U - BNMASS (4) )
                  QNORM  = DEUNEU + PRNENE
*  |  |  |  |  If we cannot split then return
                  IF ( QNORM .LE. ZERZER ) RETURN
                  CALL GRNDM(RNDM,1)
                  V = RNDM (1)
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split or into a deuteron and a neutron
*  |  |  |  |  | or into two protons and one neutron: no account is
*  |  |  |  |  | made for Coulomb effects, probability is simply assumed
*  |  |  |  |  | proportional to reaction Qs
                  IF ( V .LT. DEUNEU / QNORM ) THEN
*  |  |  |  |  | A deuteron and a proton selected
                     JEMISS = 1
                     JRESID = 3
                     RNMASS = ZMASS (3)
                     U = ZERZER
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split into 1 proton and 2 neutrons: part of the exci-
*  |  |  |  |  | tation energy is conserved to allow the further
*  |  |  |  |  | splitting of the deuteron
                  ELSE
                     JEMISS = 1
                     JRESID = 0
                     FACT = ONEONE
*  |  |  |  |  |  +----------------------------------------------------*
*  |  |  |  |  |  | Loop to compute the residual excitation energy
 1550                CONTINUE
                        FACT = FACT * 0.6666666666666667D+00
*  |  |  |  |  |  | Erncm, Eepcm are the total energies of the residual
*  |  |  |  |  |  | nucleus and of the emitted particle in the CMS frame
                        U      = FACT * PRNENE + BNMASS (3)
                        RNMASS = ZMASS (3) + U
                        ERNCM  = HLFHLF * ( UMO2 + RNMASS**2
     &                         - Z2MASS (JEMISS) ) / UMO
                        EEPCM  = UMO - ERNCM
                     IF ( EEPCM .LE. ZMASS (JEMISS) ) GO TO 1550
*  |  |  |  |  |  |
*  |  |  |  |  |  +----------------------------------------------------*
                     GO TO 2600
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |  Deuteron:
 1600          CONTINUE
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Split into a proton and a neutron if it is possible
                  IF ( U .GT. BNMASS (3) ) THEN
                     JEMISS = 1
                     JRESID = 2
                     RNMASS = ZMASS (2)
                     U = ZERZER
                     GO TO 2500
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  | Energy too low to split the deuteron, return
                  ELSE
                     WRITE (LUNERR,*)' **Dres: energy too low to split',
     &                               ' a deuteron! M2,M3',M2,M3
                     RETURN
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 1700    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
 2000    CONTINUE
         A = JA
         Z = JZ
         Q (0)  = ZERZER
         ENERG0 = FKENER (A,Z)
*  |  +----------------------------------------------------------------*
*  |  |   Note that Q(i) are not the reaction Qs but the remaining
*  |  |   energy after the reaction
         DO 2100 K = 1, 6
            JJA = JA - IA (K)
            JJZ = JZ - IZ (K)
            JJN = JJA - JJZ
            IF ( JJN .LT. 0 .OR. JJZ .LT. 0 ) THEN
               Q (K) = Q (K-1)
               GO TO 2100
            END IF
            DDJJA = JJA
            DDJJZ = JJZ
            Q (K) = MAX ( U + ENERG0 - FKENER ( DDJJA, DDJJZ )
     &            - EXMASS (K), ZERZER ) + Q (K-1)
 2100    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
*  |  +----------------------------------------------------------------*
*  |  |  If no emission channel is open then return
         IF ( Q (6) .LE. ZERZER ) THEN
            HEVSUM = SMOM1(3) + SMOM1(5) + SMOM1(6) + SMOM1(4)
            RETURN
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         CALL GRNDM(RNDM,1)
         V = RNDM (1)
         FACT = ONEONE
*  |  +----------------------------------------------------------------*
*  |  |
         DO 2200 J = 1, 6
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            IF ( V .LT. Q (J) / Q (6) ) THEN
               JEMISS = J
               JJA    = JA - IA (JEMISS)
               JJZ    = JZ - IZ (JEMISS)
*  |  |  |  +----------------------------------------------------------*
*  |  |  |  |
               DO 2150 JJ = 1, 6
*  |  |  |  |  +-------------------------------------------------------*
*  |  |  |  |  |
                  IF ( JJA .EQ. IA (JJ) .AND. JJZ .EQ. IZ (JJ) ) THEN
                     JRESID = JJ
                     RNMASS = ZMASS (JRESID)
                     ERNCM  = HLFHLF * ( UMO2 + Z2MASS (JRESID)
     &                      - Z2MASS (JEMISS) ) / UMO
                     EEPCM  = UMO - ERNCM
                     U = ZERZER
                     GO TO 2600
                  END IF
*  |  |  |  |  |
*  |  |  |  |  +-------------------------------------------------------*
 2150          CONTINUE
*  |  |  |  |
*  |  |  |  +----------------------------------------------------------*
               AJJA   = JJA
               ZJJZ   = JJZ
               RNMAS0 = AJJA * AMUMEV + FKENER ( AJJA, ZJJZ )
               GO TO 2300
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
 2200    CONTINUE
*  |  |
*  |  +----------------------------------------------------------------*
         WRITE ( LUNOUT,* )' **** error in Dres, few nucleon treatment',
     &                     ' ****'
         WRITE ( LUNERR,* )' **** error in Dres, few nucleon treatment',
     &                     ' ****'
         RETURN
*  |  +----------------------------------------------------------------*
*  |  | Loop to compute the residual excitation energy
 2300    CONTINUE
            FACT = FACT * AJJA / A
            U = FACT * ( Q (JEMISS) - Q (JEMISS-1) )
*  |  | Erncm, Eepcm are the total energies of the residual
*  |  | nucleus and of the emitted particle in the CMS frame
            RNMASS = RNMAS0 + U
            ERNCM  = HLFHLF * ( UMO2 + RNMASS**2
     &             - Z2MASS (JEMISS) ) / UMO
            EEPCM  = UMO - ERNCM
         IF ( EEPCM .LE. ZMASS (JEMISS) ) THEN
            IF ( Q (JEMISS) - Q (JEMISS-1) .GE. 1.D-06 ) GO TO 2300
*  |  +--<--<--<--<--<--< Loop back
*  |  |  Actually there is no excitation energy available!
            U = ANGLGB
            RNMASS = ONEPLS * RNMAS0
            ERNCM  = ONEPLS * RNMASS
            EEPCM  = ONEPLS * ZMASS (JEMISS)
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         GO TO 2600
*  |  From here standard two bodies kinematics with Jemiss, Rnmass
*  |  (new excitation energy is U)
 2500    CONTINUE
*  |  Erncm, Eepcm are the total energies of the residual
*  |  nucleus and of the emitted particle in the CMS frame
         ERNCM = HLFHLF * ( UMO2 + RNMASS**2 - Z2MASS (JEMISS) ) / UMO
         EEPCM = UMO - ERNCM
 2600    CONTINUE
*  |  C(i) are the direction cosines of the emitted particle
*  |  (Jemiss) in the CMS frame, of course - C(i)
*  |  are the ones of the residual nucleus (Jresid if one of the
*  |  standard partcles, say the proton)
         CALL RACO (C(1),C(2),C(3))
         PCMS  = SQRT ( EEPCM**2 - Z2MASS (JEMISS) )
*  |  Now we perform the Lorentz transformation back to the original
*  |  frame (lab frame)
*  |  First the emitted particle:
         ETAX = ETACM * COSLBR (1)
         ETAY = ETACM * COSLBR (2)
         ETAZ = ETACM * COSLBR (3)
         PCMSX = PCMS * C (1)
         PCMSY = PCMS * C (2)
         PCMSZ = PCMS * C (3)
         ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
         EPS = GAMCM * EEPCM + ETAPCM - ZMASS (JEMISS)
         PHELP = ETAPCM / ( GAMCM + ONEONE ) + EEPCM
         PLBPX = PCMSX + ETAX * PHELP
         PLBPY = PCMSY + ETAY * PHELP
         PLBPZ = PCMSZ + ETAZ * PHELP
         PHELP = SQRT (PLBPX * PLBPX + PLBPY * PLBPY + PLBPZ * PLBPZ)
         COSLBP (1) = PLBPX / PHELP
         COSLBP (2) = PLBPY / PHELP
         COSLBP (3) = PLBPZ / PHELP
*  |  Then the residual nucleus ( for it c (i) --> - c (i) ):
         EREC  = GAMCM * ERNCM - ETAPCM - RNMASS
         EKRES = 1.D-03 * EREC
         PHELP = - ETAPCM / ( GAMCM + ONEONE ) + ERNCM
         PXRES = 1.D-03 * ( - PCMSX + ETAX * PHELP )
         PYRES = 1.D-03 * ( - PCMSY + ETAY * PHELP )
         PZRES = 1.D-03 * ( - PCMSZ + ETAZ * PHELP )
         P2RES = PXRES * PXRES + PYRES * PYRES + PZRES * PZRES
         PTRES = SQRT (P2RES)
         COSLBR (1) = PXRES / PTRES
         COSLBR (2) = PYRES / PTRES
         COSLBR (3) = PZRES / PTRES
*  |  Score the emitted particle
         NPART (JEMISS) = NPART (JEMISS) + 1
         SMOM1 (JEMISS) = SMOM1 (JEMISS) + EPS
         ITEMP=NPART(JEMISS)
         EPART(ITEMP,JEMISS)=EPS
         COSEVP(1,ITEMP,JEMISS)=COSLBP(1)
         COSEVP(2,ITEMP,JEMISS)=COSLBP(2)
         COSEVP(3,ITEMP,JEMISS)=COSLBP(3)
*  |  +----------------------------------------------------------------*
*  |  |  Check if the residual nucleus is one of the emitted particles
         IF ( JRESID .GT. 0 ) THEN
            J = JRESID
            GO TO 6102
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
         JA = JA - IA (JEMISS)
         JZ = JZ - IZ (JEMISS)
* Umo is the invariant mass of the system!!
         UMO  = RNMASS
         UMO2 = UMO * UMO
         ELBTOT  = RNMASS + EREC
         GAMCM   = ELBTOT / RNMASS
         ETACM   = 1.D+03 * PTRES / RNMASS
         GO TO 1000
      END IF
*  |
*  +-------------------------------------------------------------------*
* Come to 23 at the beginning and after the end of a "normal"
* evaporation cycle
   23 CONTINUE
      A = JA
      Z = JZ
      IF(JA-8)8486,8488,8486
 8488 CONTINUE
      IF(JZ-4)8486,1224,8486
 8486 CONTINUE
      DO 2 K=1,6
         IF((A-FLA(K)).LE.(Z-FLZ(K))) THEN
            Q(K)=99999.D0
            GO TO 2
         END IF
         IF(A-TWOTWO*FLA(K))2,727,727
  727    CONTINUE
         IF(Z-TWOTWO*FLZ(K))2,728,728
  728    CONTINUE
         Q(K) = QNRG(A-FLA(K),Z-FLZ(K),A,Z) + EXMASS(K)
C 728    Q(K)=FKENER(A-FLA(K),Z-FLZ(K))-FKENER(A,Z)+EXMASS(K)
    2 CONTINUE
      FLKCOU(1)=ZERZER
      FLKCOU(2)=DOST(1,Z-FLZ(2))
      FLKCOU(3)=FLKCOU(2)+.06D+00
      FLKCOU(4)=FLKCOU(2)+.12D+00
      FLKCOU(6)=DOST(2,Z-FLZ(6))
      FLKCOU(5)=FLKCOU(6)-.06D+00
      CCOUL(1)=ONEONE
      CCOU2=DOST(3,Z-FLZ(2))
      CCOUL(2)=CCOU2+ONEONE
      CCOUL(3)=CCOU2*1.5D0+THRTHR
      CCOUL(4)=CCOU2+THRTHR
      CCOUL(6)=DOST(4,Z-FLZ(6))*TWOTWO+TWOTWO
      CCOUL(5)=TWOTWO*CCOUL(6)-ONEONE
      SIGMA=ZERZER
*  Initialize the flag which checks for open particle decay with
*  zero excitation and pairing --> for particle unstable residual
*  nuclei
      LOPPAR = .FALSE.
      DO 33 J=1,6
         IF(A-TWOTWO*FLA(J))5,725,725
  725    CONTINUE
         IF(Z-TWOTWO*FLZ(J))5,726,726
  726    CONTINUE
         MM=JA-IA(J)
         ZZ=Z-FLZ(J)
         AA=A-FLA(J)
         IF(AA.LE.ZZ)GO TO 5
*  Energy threshold for the emission of the jth-particle
         THRESH(J)=Q(J)+.88235D+00*FLKCOU(J)*FLZ(J)*
     &   ZZ/(RMASS(MM)+RHO(J))
         LOPPAR = LOPPAR .OR. THRESH (J) .GT. ZERZER
         IAA = NINT(AA)
         IZZ = NINT(ZZ)
         NN  = IAA - IZZ
*  The residual nucleus excitation energy ranges from 0 up
*  to U - Q (J)
         ILVMOD = IB0
         UMXRES = U - THRESH (J)
*  This is the a lower case of the level density
         SMALLA (J) = GETA ( UMXRES, IZZ, NN, ILVMOD, ISDUM, ASMMAX,
     &                       ASMMIN )
         CALL EEXLVL (IAA,IZZ,EEX1ST,EEX2ND,CORR)
         EEX1ST = 1.D+03 * EEX1ST
         EEX2ND = 1.D+03 * EEX2ND
         CORR   = 1.D+03 * MAX ( CORR, ZERZER )
         IF ( NN .EQ. 4 .AND. IZZ .EQ. 4 ) THEN
            IF ( U - THRESH (J) - 6.1D+00 .GT. ZERZER ) THEN
               CORR = SIXSIX
            ELSE
               TMPVAR = U-THRESH(J)-0.1D+00
               CORR = MAX ( ZERZER, TMPVAR )
            END IF
         END IF
         IF (NINT(FKEY).EQ.1) CORR=ZERZER
         CORRRR(J)=CORR
* Standard calculation:
         ARG=U-THRESH(J)-CORR
         IF(ARG)5,4,4
    5    CONTINUE
         R(J)=ZERZER
         S(J)=ZERZER
         SOS(J)=ZERZER
         GO TO 33
    4    CONTINUE
         S(J)=SQRT (SMALLA(J)*ARG)*TWOTWO
         SOS(J)=TENTEN*S(J)
   33 CONTINUE
      N1=1
      SES = MAX (S(1),S(2),S(3),S(4),S(5),S(6))
      DO 1131 J=1,6
         JS  = SOS(J) + ONEONE
         FJS = JS
         STRUN(J) = FJS - ONEONE
         IF ( S(J) .GT. ZERZER ) THEN
            MM  = JA-IA(J)
            EXPSAS=MIN(EXPMAX,MAX(EXPMIN,S(J)-SES))
            SAS = EXP (EXPSAS)
            EXPSUS=MIN(EXPMAX,MAX(EXPMIN,-S(J)))
            SUS = EXP (EXPSUS)
            EYE1(J) = ( TWOTWO * S(J)**2 -SIXSIX * S(J)
     &              + SIXSIX + SUS * ( S(J)**2 - SIXSIX ) )
     &              / ( EIGEIG * SMALLA(J)**2 )
            IF ( J .EQ. 1 ) THEN
               EYE0(J) = ( S(J) - ONEONE + SUS ) / ( TWOTWO*SMALLA(J) )
* Standard calculation
               R   (J) = RMASS(MM)**2 * ALPH(MM) * ( EYE1(J) + BET(MM)
     &                 * EYE0(J) ) * SAS
            ELSE
               R   (J) = CCOUL(J) * RMASS(MM)**2 * EYE1(J) * SAS
            END IF
            R (J) = MAX ( ZERZER, R (J) )
            SIGMA = SIGMA + R (J)
         END IF
 1131 CONTINUE
      NCOUNT = 0
 6202 CONTINUE
      IF(SIGMA)9,9,10
 9    CONTINUE
      DO 6100 J = 1,6
         IF(JA-IA(J))6100,6101,6100
 6101    CONTINUE
         IF(JZ-IZ(J))6100,6102,6100
 6100 CONTINUE
      GO TO 6099
 6102 CONTINUE
      IF ( U .GT. ANGLGB ) GO TO 1000
      JEMISS = J
C*****STORE,RESIDUAL NUC IS OF EMITTED PARTICLE TYPE
* If we are here this means that the residual nucleus is equal to
* one of the six emitted particle (the j-th one). So give to it
* all the energy, score it and return with 0 recoil and excitation
* energy for the residual nucleus
      EPS = EREC
      NPART(JEMISS) = NPART(JEMISS)+1
      ITEMP=NPART(JEMISS)
      NRNEEP = NRNEEP + 1
      SMOM1(JEMISS) = SMOM1(JEMISS) + EPS
      ITEMP=NPART(JEMISS)
      EPART(ITEMP,JEMISS)=EPS
      COSEVP(1,ITEMP,JEMISS) = COSLBR(1)
      COSEVP(2,ITEMP,JEMISS) = COSLBR(2)
      COSEVP(3,ITEMP,JEMISS) = COSLBR(3)
      GO TO 8002
*
*-->-->-->-->--> go to return
 6099 CONTINUE
      IF(JA-8)72,51,72
   51 CONTINUE
      IF(JZ-4)72,1224,72
* Come here for a "normal" step
   10 CONTINUE
      LOPPAR = .FALSE.
      CALL GRNDM(RNDM,1)
      URAN=RNDM(1)*SIGMA
      SUM = ZERZER
      DO 16 J=1,6
         K   = J
         SUM = R(J)+SUM
         IF ( SUM - URAN .GT. 0.D+00 ) GO TO 17
   16 CONTINUE
   17 CONTINUE
      JEMISS=K
      NPART(JEMISS) = NPART (JEMISS) + 1
      JS = SOS (JEMISS) + ONEONE
      IF(JS-1000)4344,4345,4345
 4345 CONTINUE
      RATIO2=(S(JEMISS)**3-6.D0*S(JEMISS)**2+15.D0*
     &S(JEMISS)-15.D0)/((2.D0*S(JEMISS)**2-6.D0*S(JEMISS)+6.D0)*SMALLA
     &(JEMISS))
      GO TO 4347
 4344 CONTINUE
      RATIO2=(P2(JS)+(P2(JS+1)-P2(JS))*
     &(SOS(JEMISS)-STRUN(JEMISS)))/SMALLA(JEMISS)
 4347 CONTINUE
      EPSAV=RATIO2*TWOTWO
*  +-------------------------------------------------------------------*
*  |  Neutron channel selected:
      IF (JEMISS .EQ. 1) THEN
         MM=JA-IA(J)
         EPSAV=(EPSAV+BET(MM))/(ONEONE+BET(MM)*EYE0(JEMISS)
     &        /EYE1(JEMISS))
*  |  +----------------------------------------------------------------*
*  |  |  Compute the fission width relative to the neutron one:
*  |  |  this part is taken from subroutine EMIT of LAHET
         IF ( IFISS .GT. 0 .AND. JZ .GE. 71 .AND. .NOT. FISINH ) THEN
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Compute the correction factor for the fission width:
            IF ( JZ .GT. 88 ) THEN
               AGOES = ONEONE
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  |
            ELSE
               AGOES = MAX ( ONEONE, ( U-SEVSEV ) / ( EPSAV+SEVSEV ) )
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
*  |  |  Finally this is the relative fission width:
*  |  |  This is : Probfs = 1 / ( 1 + G_n / G_f )
            PROBFS = FPROB ( Z, A, U ) / AGOES
*  |  |  +-------------------------------------------------------------*
*  |  |  |  Check if it will be fission:
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. PROBFS ) THEN
               FISINH = .TRUE.
               KFISS  = 1
*  |  |  |  Undo the counting of the neutron evaporation
               NPART (JEMISS) = NPART (JEMISS) -1
               GO TO 9260
            END IF
*  |  |  |
*  |  |  +-------------------------------------------------------------*
         END IF
*  |  |
*  |  +----------------------------------------------------------------*
      END IF
*  |
*  +-------------------------------------------------------------------*
  20  CONTINUE
      CALL GRNDM(RNDM,2)
      E1=-HLFHLF*LOG(RNDM(1))
      E2=-HLFHLF*LOG(RNDM(2))
* Eps should be the total kinetic energy in the CMS frame
* Standard calculation:
      EPS=(E1+E2)*EPSAV+THRESH(JEMISS)-Q(JEMISS)
      AR = A - IA(JEMISS)
      ZR = Z - IZ(JEMISS)
* The CMS energy is updated
      IMASS = NINT (AR)
      IF ( IMASS .EQ. 8 .AND. NINT (ZR) .EQ. 4 ) THEN
         UNEW = U - EPS - Q(JEMISS)
         UMAX = U - THRESH(JEMISS)
         IF ( UNEW .GT. 6.D+00 ) THEN
            UMIN = 6.D+00
         ELSE IF ( UNEW .GT. 4.47D+00 .AND. UMAX .GT. 6.D+00 ) THEN
            UMIN = 4.47D+00
            UNEW = 6.D+00
         ELSE IF ( UNEW .GT. 1.47D+00 .AND. UMAX .GT. 2.94D+00 ) THEN
            UMIN = 1.47D+00
            UNEW = 2.94D+00
         ELSE
            UMIN = -0.1D+00
            UNEW = ANGLGB * 0.1D+00
         END IF
      ELSE IF ( IMASS .LE. 4 ) THEN
         IPRO = NINT ( ZR )
         INEU = IMASS - IPRO
         IF ( IMASS .EQ. 1 ) THEN
*  Be sure that residual neutrons or protons are not left excited
            UMIN = 0.D+00
            UNEW = 0.D+00
            EPS  = U - Q(JEMISS)
         ELSE IF ( IPRO .EQ. 0 .OR. INEU .EQ. 0 ) THEN
*  Ipro protons or ineu neutrons arrived here!
            UMIN = CORRRR(JEMISS)
            UNEW = U - EPS - Q(JEMISS)
         ELSE IF ( IMASS .LE. 2 ) THEN
*  Be sure that residual deuterons are not left excited!
            UMIN = 0.D+00
            UNEW = 0.D+00
            EPS  = U - Q(JEMISS)
         ELSE IF ( ABS ( INEU - IPRO ) .LE. 1 ) THEN
*  For the moment also residual 3-H, 3-He and 4-He are not left
*  excited !
            UMIN = 0.D+00
            UNEW = 0.D+00
            EPS  = U - Q(JEMISS)
         ELSE
            UMIN = CORRRR(JEMISS)
            UNEW = U - EPS - Q(JEMISS)
         END IF
      ELSE
         UMIN = CORRRR(JEMISS)
         UNEW = U - EPS - Q(JEMISS)
      END IF
* Standard calculation
      IF(UNEW-UMIN)6200,6220,6220
 6220 CONTINUE
*or   RNMASS = AR * AMUMEV + FKENER (AR,ZR)
      RNMASS = AR * AMUMEV + FKENER (AR,ZR) + UNEW
      UMIN2  = ( RNMASS + ZMASS (JEMISS) )**2
      IF ( UMIN2 .GE. UMO2 ) THEN
         GO TO 6200
      END IF
      U = UNEW
* C(i) are the direction cosines of the evaporated particle in the CMS
* frame, of course - C(i) are the ones of the residual nucleus
      CALL RACO(C(1),C(2),C(3))
* Erncm, Eepcm are the total energies of the residual nucleus and
* of the evaporated particle in the CMS frame
      ERNCM = 0.5D+00 * ( UMO2 + RNMASS**2 - Z2MASS (JEMISS) ) / UMO
      EEPCM = UMO - ERNCM
      PCMS  = SQRT ( EEPCM**2 - Z2MASS (JEMISS) )
* Now we perform the Lorentz transformation back to the original
* frame (lab frame)
* First the evaporated particle:
      ETAX = ETACM * COSLBR (1)
      ETAY = ETACM * COSLBR (2)
      ETAZ = ETACM * COSLBR (3)
      PCMSX = PCMS * C (1)
      PCMSY = PCMS * C (2)
      PCMSZ = PCMS * C (3)
      ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
      EPS = GAMCM * EEPCM + ETAPCM - ZMASS (JEMISS)
      PHELP = ETAPCM / (GAMCM + 1.D0) + EEPCM
      PLBPX = PCMSX + ETAX * PHELP
      PLBPY = PCMSY + ETAY * PHELP
      PLBPZ = PCMSZ + ETAZ * PHELP
      PHELP = SQRT (PLBPX * PLBPX + PLBPY * PLBPY + PLBPZ * PLBPZ)
      COSLBP (1) = PLBPX / PHELP
      COSLBP (2) = PLBPY / PHELP
      COSLBP (3) = PLBPZ / PHELP
* Then the residual nucleus ( for it c (i) --> - c (i) ):
      EREC  = GAMCM * ERNCM - ETAPCM - RNMASS
      EKRES = 1.D-03 * EREC
      PHELP = - ETAPCM / (GAMCM + 1.D0) + ERNCM
      PXRES = 1.D-03 * ( - PCMSX + ETAX * PHELP )
      PYRES = 1.D-03 * ( - PCMSY + ETAY * PHELP )
      PZRES = 1.D-03 * ( - PCMSZ + ETAZ * PHELP )
      P2RES = PXRES * PXRES + PYRES * PYRES + PZRES * PZRES
      PTRES = SQRT (P2RES)
      COSLBR (1) = PXRES / PTRES
      COSLBR (2) = PYRES / PTRES
      COSLBR (3) = PZRES / PTRES
* Check energy and momentum conservation !!
      IF (EREC .LE. 0.D+00) THEN
         PTRES = 0.D+00
         EREC  = 0.D+00
      END IF
* Umo is the invariant mass of the system!!
      UMO  = RNMASS
      UMO2 = UMO * UMO
      ELBTOT  = RNMASS + EREC
      GAMCM   = ELBTOT / RNMASS
      ETACM   = 1.D+03 * PTRES / RNMASS
      GO TO 76
 6200 CONTINUE
      NCOUNT = NCOUNT + 1
      IF ( NCOUNT .LE. 10 ) GO TO 20
      SIGMA = SIGMA - R(JEMISS)
* if we are here we have sampled for > 10 times a negative energy Unew
      NPART(JEMISS)=NPART(JEMISS)-1
      R(JEMISS) = 0.D0
      NCOUNT = 0
      GO TO 6202
   76 CONTINUE
      JAT=JA-IA(JEMISS)
      JZT=JZ-IZ(JEMISS)
      IF(JAT-JZT)9,9,77
   77 CONTINUE
      JA=JAT
      JZ=JZT
C*****STORE,END OF NORMAL CYCLE
      SMOM1(JEMISS)=SMOM1(JEMISS)+EPS
      ITEMP=NPART(JEMISS)
      EPART(ITEMP,JEMISS)=EPS
      COSEVP(1,ITEMP,JEMISS)=COSLBP(1)
      COSEVP(2,ITEMP,JEMISS)=COSLBP(2)
      COSEVP(3,ITEMP,JEMISS)=COSLBP(3)
* The following card switch to the rough splitting treatment
      IF (JA .LE. 2) GO TO 1000
      IF(JA-8)23,1223,23
 1223 CONTINUE
      IF(JZ-4)23,1224,23
* If we are here the residual nucleus is a 8-Be one, break it into
* two alphas with all the available energy (U plus the Q of the breakup)
* , score them and return with 0 recoil and excitation energy
 1224 CONTINUE
      LOPPAR = .FALSE.
      IF(U)1228,1229,1229
 1228 CONTINUE
      EPS=0.D0
      GO TO 1230
 1229 CONTINUE
 1230 CONTINUE
      NBE8=NBE8+1
* C(i) are the direction cosines of the first alpha in the CMS
* frame, of course - C(i) are the ones of the other
      CALL RACO(C(1),C(2),C(3))
* Ecms is the total energy of the alphas in the CMS frame
      ECMS  = 0.5D+00 * UMO
      PCMS  = SQRT ( ECMS**2 - Z2MASS (6) )
* Now we perform the Lorentz transformation back to the original
* frame (lab frame)
* First alpha:
      ETAX = ETACM * COSLBR (1)
      ETAY = ETACM * COSLBR (2)
      ETAZ = ETACM * COSLBR (3)
      PCMSX = PCMS * C (1)
      PCMSY = PCMS * C (2)
      PCMSZ = PCMS * C (3)
      ETAPCM = PCMSX * ETAX + PCMSY * ETAY + PCMSZ * ETAZ
      EPS = GAMCM * ECMS + ETAPCM - ZMASS (6)
      PHELP = ETAPCM / (GAMCM + 1.D0) + ECMS
      PLBPX = PCMSX + ETAX * PHELP
      PLBPY = PCMSY + ETAY * PHELP
      PLBPZ = PCMSZ + ETAZ * PHELP
      PHELP = SQRT (PLBPX * PLBPX + PLBPY * PLBPY + PLBPZ * PLBPZ)
* Store the first alpha!!
      SMOM1(6) = SMOM1(6) + EPS
      NPART(6) = NPART(6) + 1
      ITEMP = NPART(6)
      EPART (ITEMP,6) = EPS
      COSEVP (1,ITEMP,6) = PLBPX / PHELP
      COSEVP (2,ITEMP,6) = PLBPY / PHELP
      COSEVP (3,ITEMP,6) = PLBPZ / PHELP
* Then the second alpha ( for it c (i) --> - c (i) ):
      EPS = GAMCM * ECMS - ETAPCM - ZMASS (6)
      PHELP = - ETAPCM / (GAMCM + 1.D0) + ECMS
      PLBPX = - PCMSX + ETAX * PHELP
      PLBPY = - PCMSY + ETAY * PHELP
      PLBPZ = - PCMSZ + ETAZ * PHELP
      PHELP = SQRT (PLBPX * PLBPX + PLBPY * PLBPY + PLBPZ * PLBPZ)
* Store the second alpha !!
      SMOM1(6) = SMOM1(6) + EPS
      NPART(6) = NPART(6) + 1
      ITEMP = NPART(6)
      EPART (ITEMP,6) = EPS
      COSEVP (1,ITEMP,6) = PLBPX / PHELP
      COSEVP (2,ITEMP,6) = PLBPY / PHELP
      COSEVP (3,ITEMP,6) = PLBPZ / PHELP
 8002 CONTINUE
      LOPPAR = .FALSE.
      EREC   = ZERZER
      U      = ZERZER
      EKRES  = ZERZER
      PTRES  = ZERZER
   72 CONTINUE
      HEVSUM=SMOM1(3)+SMOM1(5)+SMOM1(6)+SMOM1(4)
      RETURN
*.......................................................................
*///// RAL FISSION ROUTINE /////
 9260 CONTINUE
*  +-------------------------------------------------------------------*
*  |  Record the direction cosines of the fissioning nucleus
      DO 9270 I=1,3
         COSLF0 (I) = COSLBR (I)
 9270 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      CALL FISFRA ( JA, JZ, U, EREC, UMO, GAMCM, ETACM )
*  +-------------------------------------------------------------------*
*  |  Check for fission failures!!
      IF ( .NOT. FISINH ) THEN
         PENBAR = .FALSE.
         GO TO 23
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Do not pick up the fission fragments, rather go back to Evevap
      HEVSUM=SMOM1(3)+SMOM1(5)+SMOM1(6)+SMOM1(4)
      RETURN
      END
