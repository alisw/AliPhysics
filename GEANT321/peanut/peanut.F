*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:02  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE PEANUT.FOR
*COPY PEANUT
*
*=== peanut ===========================================================*
*
      SUBROUTINE PEANUT ( KPROJ, EKE, PPROJ, TXX, TYY, TZZ, WEE )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
*
#include "geant321/balanc.inc"
#include "geant321/eva0.inc"
#include "geant321/fheavy.inc"
#include "geant321/finuc.inc"
#include "geant321/higfis.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
*
      COMMON / FKCOSP / C1ST (3), C2ND (3), LEMISS
      LOGICAL LEMISS
      COMMON / FKCMCY / NPCYCL (20,2), IEVT, LOUT
*
      COMMON / FKPLOC / IABCOU
      LOGICAL LBCHCK, LBIMPC, LTRPPD, LASSOR, LEXIT, LPRCYC, LEXPLC,
     &        LNWINT
      DIMENSION IPTYPE (39)
      REAL RNDM(2)
      SAVE IPTYPE
      DATA IPTYPE /  1,  2,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,
     &               3,  3,  4,  4,  5,  5,  0,  5,  5,  5,  3,  4,
     &               4,  0,  0,  0,  0,  0,  6,  6,  6,  7,  8,  7,
     &               8,  9, 10 /
*
*
      IEVPRE = IEVPRE + 1
      NUSCIN = 0
      IABCOU = 0
      IF ( EKE .GT. 2.D+00 * GAMMIN ) THEN
         EOTEST = ETTOT
         PTORI  = PPROJ
         PXORI  = PTORI * TXX
         PYORI  = PTORI * TYY
         PZORI  = PTORI * TZZ
      ELSE
         EOTEST = ETTOT - EKE
         ETTOT  = EOTEST
         PTORI  = 0.D+00
         PXTTOT = 0.D+00
         PYTTOT = 0.D+00
         PZTTOT = 0.D+00
         PTTOT  = 0.D+00
         PXORI  = 0.D+00
         PYORI  = 0.D+00
         PZORI  = 0.D+00
      END IF
      ETEPS  = 1.D-10 * ETTOT
      ICHTOT = ICHTAR + ICH  (KPROJ)
      IBTOT  = IBTAR  + IBAR (KPROJ)
      IBNUCL = IBTAR
      IBORI  = IBAR   (KPROJ)
      IPTORI = IPTYPE (KPTOIP(KPROJ))
      KPORI  = KPROJ
      EKORI  = EKE
      ZZTAR  = ICHTAR
      BBTAR  = IBTAR
      IF ( ICH (KPROJ) .NE. 0 .AND. EKE .GT. 2.D+00 * GAMMIN ) THEN
         FLKCOU = DOST ( 1, ZZTAR )
         CCOUL  = DOST ( 3, ZZTAR )
         CLMBBR = ICH (KPROJ) * COULBH * ZZTAR / RMASS (IBTAR)
         IF ( CLMBBR .GT. 0.9D+00 * EKE ) THEN
            TMPEKE = 0.98D+00 * EKE
            CLMHLP = MIN ( CLMBBR * FLKCOU, TMPEKE )
            CLMBBR = MIN  ( CLMBBR, EKE )
            WEIGH1 = 10.0D+00 * ( CLMBBR / EKE - 0.9D+00 )
            CLMBBR = WEIGH1 * CLMHLP + ( 1.D+00 - WEIGH1 ) * CLMBBR
         END IF
         BFCLMB = SQRT ( 1.D+00 - CLMBBR / EKE )
         RDCLMB = ICH (KPROJ) * COULPR * ZZTAR / CLMBBR
      ELSE
         CLMBBR = 0.D+00
         BFCLMB = 1.D+00
         RDCLMB = AINFNT
      END IF
      IBRES = IBTOT
      ICRES = ICHTOT
      BBRES = IBRES
      ZZRES = ICRES
      AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER ( BBRES, ZZRES )
      AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
      NPROT = 0
      NNEUT = 0
      NHOLE = 0
      AVEBIN = ( ( BBTAR - ZZTAR ) * AMNUCL (2) + ZZTAR * AMNUCL (1)
     &       - AMNTAR ) / BBTAR
      AMMHLP = ( BBTAR - 1.D+00 ) * AMUAMU + 1.D-03
     &       * FKENER ( BBTAR - ONEONE, ZZTAR - ONEONE )
      AMNHLP = AMMHLP - ( ZZTAR - 1.D+00 ) * AMELEC + ELBNDE (ICHTAR-1)
      BNENRG (1) = AMNHLP + AMNUCL (1) - AMNTAR
      AMMHLP = ( BBTAR - 1.D+00 ) * AMUAMU + 1.D-03
     &       * FKENER ( BBTAR - ONEONE, ZZTAR )
      AMNHLP = AMMHLP - ZZTAR * AMELEC + ELBNDE (ICHTAR)
      BNENRG (2) = AMNHLP + AMNUCL (2) - AMNTAR
 
      BNENRG (3) = 0.5D+00 * ( BNENRG (1) + BNENRG (2) )
      RHORED = 1.D+00
      NPNUC  = 0
      NNUCTS = 0
      NHLEXP = 0
      JNUCTS = 0
      ACOLL  = ANOW
      ZCOLL  = ZNOW
      IF ( .NOT. LPREEX ) THEN
         IF ( KPROJ .EQ. 1 .OR. KPROJ .EQ. 8 ) THEN
            LEXPLC = EKE .GT. 0.250D+00
         ELSE IF ( KPROJ .EQ. 14 .AND. PTTOT .LE. 0.D+00 ) THEN
            LEXPLC = .TRUE.
         ELSE
            STOP 'LEXPLC'
         END IF
      ELSE
         LEXPLC = .TRUE.
      END IF
      IF ( LEXPLC .AND. EKE .GT. IBAR(KPROJ) * EKEEXP ) THEN
         ICYCL  = 0
         IREINT = 0
         LPRCYC = .TRUE.
         LBCHCK = .FALSE.
         LBIMPC = .TRUE.
         LELSTC = .FALSE.
         LABRST = .FALSE.
         LABSRP = .FALSE.
         LINELS = .FALSE.
         LCHEXC = .FALSE.
         RHOEXP = 0.D+00
         EKFEXP = 0.D+00
         EKFREI = 0.D+00
         RHOREI = 0.D+00
         KPRIN  = KPROJ
         KRFLIN = 0
         ERECLD = 0.D+00
         BNPREV = 0.D+00
         EKECON = EKE
         PNUCCO = PPROJ
         CALL BIMSEL ( KPROJ, TXX, TYY, TZZ, LBCHCK )
         WLLPRO = WLLRED
         BNPROJ = WLLRED * BNDNUC
         RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
         EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
         IAAFT  = IBRES - IBAR (KPROJ)
         IZAFT  = ICRES - ICH  (KPROJ)
         CALL EEXLVL ( IAAFT, IZAFT, EEXDEL, EEXMIN, EEXDUM )
         EEXANY = EEXDEL
         IF ( ALPHAL .LE. 0.D+00 ) THEN
            DEFPRO = 0.D+00
            DEFNEU = 0.D+00
         ELSE
            DEFPRO = 0.D+00
            DEFNEU = 0.D+00
         END IF
         DEFNUC (1) = DEFPRO
         DEFNUC (2) = DEFNEU
  100    CONTINUE
         IF ( LELSTC ) THEN
            CALL NUCNUC ( IKPMX , KRFLIN, WEE   , ERECMN, LBIMPC,
     &                    LBCHCK, ICYCL , NHOLE , NPROT , NNEUT ,
     &                    LEXIT , LNWINT )
         ELSE IF ( KPRIN .EQ. 14 ) THEN
            CALL PIOABS ( IKPMX , KRFLIN, WEE   , ERECMN, LBIMPC,
     &                    LBCHCK, ICYCL , NHOLE , NPROT , NNEUT ,
     &                    LEXIT , LNWINT )
         ELSE
            STOP 'Int_kind'
         END IF
         IF ( LNWINT ) GO TO 100
         BBRES = IBRES
         ZZRES = ICRES
         IF ( .NOT. LEXIT ) THEN
            LPRCYC = .FALSE.
         ELSE
            BNPREV = BNPREV + BNDUSE
         END IF
  200    CONTINUE
         LELSTC = .FALSE.
         LABRST = .FALSE.
         LABSRP = .FALSE.
         LINELS = .FALSE.
         LCHEXC = .FALSE.
         GAMMAX = 0.D+00
         IREFMN = 10000
         IKPMX  = 0
         IBCHCK = 0
         ICCHCK = 0
         IBNUCL = 0
         ICNUCL = 0
         ILIVE  = 0
         LTRPPD = .FALSE.
         DO 300 KP = 1, NPNUC
            IF ( KPNUCL (KP) .LE. 0 ) GO TO 300
            ILIVE = ILIVE + 1
            KPNUC = KPNUCL (KP)
            IPTNUC = IPTYPE (KPTOIP(KPNUC))
            IF ( IPTNUC .EQ. 1 ) THEN
               BNDNU0 = BNENRG (1+KPNUC/8)
               WLLRE0 = POTBAR
            ELSE
               IF ( IBAR (KPNUC) .NE. 0 ) THEN
                  WLLRE0 = POTBAR
                  BNDNU0 = BNENRG (3)
               ELSE IF ( KPNUC .LE. 11 ) THEN
                  WLLRE0 = 0.D+00
                  BNDNU0 = 0.D+00
               ELSE
                  WLLRE0 = POTMES
                  BNDNU0 = BNENRG (3)
               END IF
            END IF
            IF ( EKFNUC (KP) .GT. -100.D+00 ) THEN
               GAMMA = ( ENNUC (KP) - WLLRE0 * ( EKFNUC (KP) + BNDNU0 )
     &               ) / AM (KPNUC)
            ELSE
               IF ( AM (KPNUC) .LE. ANGLGB ) THEN
                  GAMMA = AINFNT
               ELSE
                  GAMMA = ENNUC (KP) / AM (KPNUC)
               END IF
            END IF
            IF ( IBAR (KPNUC) .GT. 0 ) THEN
               IBNUCL = IBNUCL + IBAR (KPNUC)
               ICNUCL = ICNUCL + ICH  (KPNUC)
            END IF
            IBCHCK = IBCHCK + IBAR (KPNUC)
            ICCHCK = ICCHCK + ICH  (KPNUC)
            IF ( KRFNUC (KP) .LT. IREFMN ) THEN
               IREFMN = KRFNUC (KP)
               GAMMAX = GAMMA
               IKPMX  = KP
               WLLRED = WLLRE0
               BNDNUC = BNDNU0
            ELSE IF ( KRFNUC (KP) .EQ. IREFMN ) THEN
               IF ( GAMMA .GT. GAMMAX ) THEN
                  GAMMAX = GAMMA
                  IKPMX  = KP
                  WLLRED = WLLRE0
                  BNDNUC = BNDNU0
               END IF
            END IF
  300    CONTINUE
         IBNUCL = IBRES - IBNUCL - NPROT - NNEUT
         ICNUCL = ICRES - ICNUCL - NPROT
         ACOLL  = IBNUCL
         ZCOLL  = ICNUCL
         RHORED = ACOLL / BBTAR
         IF ( IKPMX .LE. 0 ) THEN
            IBCKC = IBTOT  - IBINTR - IBNUCR
            ICCKC = ICHTOT - ICINTR - ICNUCR
            IF ( IBCKC .NE. IBRES .OR. ICCKC .NE. ICRES ) THEN
               ICRES = ICCKC
               IBRES = IBCKC
            END IF
            NEXPEM = NP-NP0
            DO 450 IJJ = 1,IGREYN
               NPCYCL (IJJ,1) = 0
  450       CONTINUE
            DO 460 IJJ = 1,IGREYP
               NPCYCL (IJJ,2) = 0
  460       CONTINUE
            BBRES = IBRES
            ZZRES = ICRES
            ANOW  = BBRES
            ZNOW  = ZZRES
            AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER ( BBRES, ZZRES)
            AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
            PXRES = PXTTOT - PXNUCR - PXINTR
            PYRES = PYTTOT - PYNUCR - PYINTR
            PZRES = PZTTOT - PZNUCR - PZINTR
            PTRES2= PXRES**2 + PYRES**2 + PZRES**2
            PTRES = SQRT ( PTRES2 )
            ERES  = ETTOT - EINTR - ENUCR
            UMO2  = ( ERES - PTRES ) * ( ERES + PTRES )
            IF ( UMO2 .LT. ONEMNS*AMNRES**2 ) THEN
               UMO = SQRT (UMO2)
               WRITE ( LUNOUT,* )' 2:UMO,AMNRES',UMO,AMNRES
               GO TO 530
            ELSE IF ( UMO2 - AMNRES*AMNRES .LT. AMNRES*TVEPSI ) THEN
               UMO = SQRT (UMO2)
               GO TO 530
            END IF
            IF ( ICYCL .NE. NHOLE - IABCOU ) THEN
               WRITE (LUNOUT,*)' *** KPORI, ICYCL, NHOLE, IABCOU',
     &                               KPORI, ICYCL, NHOLE, IABCOU
               ICYCL = NHOLE - IABCOU
            END IF
            NPTOT  = NPROT + NNEUT
            NHLEXP = NHOLE
            IF ( .NOT. LPRCYC .OR. NPTOT .LE. 0 .OR. NNUCTS .GT. 0 )
     &         THEN
               LEMISS = .FALSE.
               NHOLE  = NHOLE + 1
               IF ( EKFREI .GT. ANGLGB ) THEN
                  RHOIMP = ( RHOEXP + RHOREI ) / NHOLE
                  EKFIMP = ( EKFEXP + EKFREI ) / NHOLE
               ELSE
                  RHOIMP = ( RHOEXP + RHOAVE ) / NHOLE
                  EKFIMP = ( EKFEXP + EKFAVE ) / NHOLE
               END IF
               ANPROT = NPROT
               ANNEUT = NNEUT
               ACOLL  = ACOLL - 1.D+00
               IF ( NPTOT .GT. 0 ) THEN
                  PNPROT = ( ZNOW - ANPROT ) * ( 3.D+00 * ANNEUT
     &                   + ANPROT ) / ( ANPROT * ( ZNOW - ANPROT
     &                   + 3.D+00 * ( ANOW - ANNEUT - ZNOW ) ) + ANNEUT
     &                   * ( 3.D+00 * ( ZNOW
     &                   - ANPROT ) + ANOW - ANNEUT - ZNOW ) )
               ELSE
                  PNPROT = ZNOW / ANOW
               END IF
               CALL GRNDM(RNDM,1)
               IF ( RNDM (1) .LT. PNPROT ) THEN
                  NPROT = NPROT + 1
                  ZCOLL = ZCOLL - 1.D+00
               ELSE
                  NNEUT = NNEUT + 1
               END IF
            ELSE
               RHOIMP = RHOEXP / NHLEXP
               EKFIMP = EKFEXP / NHLEXP
               IF ( ICYCL .EQ. 1 ) THEN
                  IF ( NP .LE. NP0 ) THEN
                     LEMISS = .FALSE.
                     ICYCL  = ICYCL - 1
                  ELSE
                     LEMISS = .TRUE.
                  END IF
               ELSE
                  LEMISS = .FALSE.
                  ICYCL  = ICYCL - 1
               END IF
            END IF
            GO TO 500
         END IF
         KPNUC  = KPNUCL (IKPMX)
         IPTNUC = IPTYPE (KPTOIP(KPNUC))
         ERECMN = MAX ( ERECLD, ERECMN )
         ERECLD = ERECMN
         ERECMN = ERECMN / ( ICYCL + IGREYP + IGREYN )
         IAAFT  = IBRES - IBAR (KPNUC)
         IZAFT  = ICRES - ICH  (KPNUC)
         CALL EEXLVL ( IAAFT, IZAFT, EEXDEL, EEXMIN, EEXDUM )
         IF ( NP .EQ. NP0 .AND. KPNUC .EQ. KPROJ ) THEN
            EEXANY = EEXDEL
         ELSE
            EEXANY = 0.D+00
         END IF
         AAFT = BBRES - IBAR (KPNUC)
         ZAFT = ZZRES - ICH  (KPNUC)
         AMMAFT = AAFT * AMUAMU + 0.001D+00 * FKENER ( AAFT, ZAFT )
         AMNAFT = AMMAFT - ZAFT * AMELEC + ELBNDE ( NINT (ZAFT) )
         IF ( WLLRED .GT. ANGLGB ) THEN
            IF ( EKFNUC (IKPMX) .GT. -100.D+00 ) THEN
               BNDGEN = IBAR (KPNUC) * AM (KPNUC) + AMNAFT - AMNRES
               IF ( NP .EQ. NP0 .AND. IPTNUC .EQ. IPTORI ) THEN
                  BNDUSE = BNPROJ + AMNAFT - AMNTAR + AM (KPNUC)
     &                   - AM (KPROJ)
               ELSE IF ( NP .EQ. NP0 ) THEN
                  IF ( IPTNUC .EQ. 1 ) THEN
                     BNDUSE = AMNAFT - AMNTAR + AM (KPNUC)
                     BNDUSE = MAX ( BNDUSE, ZERZER )
                  ELSE
                     BNDUSE = WLLRED * BNDNUC
                  END IF
               ELSE IF ( NUSCIN .EQ. 1 ) THEN
                  AMEMIT = 0.D+00
                  DO 430 KP = NP0+1, NP
                     IPTPAR = IPTYPE (KPTOIP(KPART(KP)))
                     IF ( IPTPAR .EQ. 1 .OR. IPTPAR .EQ. IPTORI )
     &               AMEMIT = AMEMIT + AM (KPTOIP(KPART(KP)))
  430             CONTINUE
                  IF ( IPTNUC .EQ. IPTORI ) THEN
                     BNTRUE = AMNAFT + AMEMIT + AM (KPNUC) - AMNTAR
     &                      - AM (KPROJ)
                     BNDUSE = BNPROJ + BNTRUE - BNPREV
                     BNDUSE = MAX ( BNDUSE, ZERZER )
                  ELSE IF ( IPTNUC .EQ. 1 ) THEN
                     BNDUSE = AMNAFT + AMEMIT + AM (KPNUC) - AMNTAR
     &                      - BNPREV
                     BNDUSE = MAX ( BNDUSE, ZERZER )
                  ELSE
                     BNDUSE = WLLRED * BNDNUC
                  END IF
               ELSE
                  BNDUSE = WLLRED * MAX ( BNDGEN, ZERZER )
               END IF
               EKFPRE = EKFNUC (IKPMX)
               VWELL0 = WLLRED * EKFNUC (IKPMX) + BNDUSE + ERECMN
               ENNUC  (IKPMX) = ENNUC (IKPMX) - VWELL0
               EKFNUC (IKPMX) = WLLRED * BNDNUC - BNDUSE - ERECMN
            ELSE
               BNDGEN = IBAR (KPNUC) * AM (KPNUC) + AMNAFT - AMNRES
               IF ( NP .EQ. NP0 .AND. KPNUC .EQ. KPORI ) THEN
                  BNDUSE = BNPROJ
               ELSE
                  BNDUSE = WLLRED * MAX ( BNDGEN, ZERZER )
               END IF
               EKFPRE = EKFAVE
               RHNUCL (IKPMX) = RHOAVE
               VWELL0 = BNDUSE - WLLRED * BNDNUC + ERECMN
               ENNUC  (IKPMX) = ENNUC (IKPMX) - VWELL0
               EKFNUC (IKPMX) = - VWELL0
            END IF
         ELSE
            VWELL0 = 0.D+00
            EKFNUC (IKPMX) = 0.D+00
            EKFPRE = 0.D+00
            BNDUSE = 0.D+00
         END IF
         EKNNUC = ENNUC (IKPMX) - AM (KPNUC)
         EKECON = EKNNUC - EKFNUC (IKPMX)
         IF ( ICH (KPNUC) .GT. 0 ) THEN
            FLKCOU = DOST ( 1, ZAFT )
            ETHRES = FLKCOU * ICH (KPNUC) * COULBH * ZAFT
     &             / RMASS ( NINT (AAFT) )
            IF ( EKNNUC .GT. ETHRES ) THEN
               LASSOR = .FALSE.
               FREJE  = 1.D+00 - ( ETHRES / EKNNUC )**3
               CALL GRNDM(RNDM,1)
               IF ( RNDM (1) .GE. FREJE ) LTRPPD = .TRUE.
            ELSE
               LASSOR = .TRUE.
            END IF
         ELSE
            ETHRES = 0.D+00
            IF ( EKNNUC .GT. ETHRES ) THEN
               LASSOR = .FALSE.
            ELSE
               LASSOR = .TRUE.
            END IF
         END IF
         IF ( LASSOR .OR. LTRPPD ) THEN
            IF ( KPNUC .EQ. 1 .OR. KPNUC .EQ. 8 ) THEN
               KPNUCL (IKPMX) = - KPNUCL (IKPMX)
               ENNUC (IKPMX) = ENNUC (IKPMX) - EKFNUC (IKPMX)
               EKFNUC (IKPMX) = -AINFNT
               NPROT  = NPROT + ICH (KPNUC)
               NNEUT  = NNEUT + 1 - ICH (KPNUC)
               IREINT = IREINT + 1
               LPRCYC = .FALSE.
               GO TO 200
            ELSE
               LTRPPD = .TRUE.
               STOP 'KPNUCL_TRAPPED'
            END IF
         ELSE
            LTRPPD = .FALSE.
         END IF
         IKPNWI = IKPMX
         IF ( KPNUC .EQ. 1 .OR. KPNUC .EQ. 8 ) THEN
            ETHMNM = ETHRES + EKEMNM
            ETHREI = MAX ( EKREXP, ETHMNM )
            IF ( EKECON .LT. ETHMNM ) THEN
               NPROT = NPROT + ICH (KPNUC)
               NNEUT = NNEUT + 1 - ICH (KPNUC)
               ENNUC  (IKPMX) = EKECON + AM (KPNUC)
               KPNUCL (IKPMX) = - KPNUCL (IKPMX)
               NNUCTS = NNUCTS + 1
               INUCTS (NNUCTS) = IKPMX
               JNUCTS = NUSCIN
               EKFNUC (IKPMX) = EKFPRE
               ENNUC  (IKPMX) = ENNUC (IKPMX) - BNDUSE + BNDGEN
               RSTNUC (IKPMX) = BNDGEN
               LPRCYC = .TRUE.
               GO TO 200
            ELSE IF ( EKECON .LT. ETHREI ) THEN
               LBCHCK = .FALSE.
               IKPNWI = - IKPMX
            ELSE
               LBCHCK = .FALSE.
            END IF
         ELSE
            LBCHCK = .FALSE.
         END IF
         PNUCCO = SQRT ( EKECON * ( EKECON + 2.D+00 * AM (KPNUC) ) )
         CALL NWISEL ( IKPNWI, LBCHCK )
  350    CONTINUE
         IF ( BIMPCT .GT. RADTOT .AND. .NOT. LTRPPD ) THEN
            KPRIN = KPNUC
            IF ( EKNNUC .NE. EKECON ) PNUCCO = SQRT ( EKNNUC * ( EKNNUC
     &         + 2.D+00 * AM (KPNUC) ) )
            IF ( ABS ( PNUCL (IKPMX) - PNUCCO ) .GT. ANGLGB * PNUCCO )
     &           CALL PHDSET ( IKPMX )
            IBRES = IBRES - IBAR (KPNUC)
            ICRES = ICRES - ICH  (KPNUC)
            BBRES = IBRES
            ZZRES = ICRES
            AMMRES = AMMAFT
            AMNRES = AMNAFT
            CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTRPPD )
            EKNNUC = ENNUC (IKPMX) - AM (KPNUC)
            IF ( LTRPPD ) GO TO 350
            NP = NP + 1
            TKI   (NP) = ENNUC  (IKPMX) - AM (KPNUC)
            KPART (NP) = KPNUC
            PLR   (NP) = PNUCL  (IKPMX)
            CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
            CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
            CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
            WEI   (NP) = WEE
            KPNUCL (IKPMX) = 0
            IF ( KPNUC .EQ. 1 .OR. KPNUC .EQ. 8 ) THEN
               IGREYP = IGREYP + ICH (KPNUC)
               IGREYN = IGREYN + 1 - ICH (KPNUC)
               PXINTR = PXINTR + PXNUCL (IKPMX)
               PYINTR = PYINTR + PYNUCL (IKPMX)
               PZINTR = PZINTR + PZNUCL (IKPMX)
               EINTR  = EINTR  + ENNUC  (IKPMX)
               IBINTR = IBINTR + IBAR   (KPNUC)
               ICINTR = ICINTR + ICH    (KPNUC)
            ELSE
               IOTHER = IOTHER + 1
               PXNUCR = PXNUCR + PXNUCL (IKPMX)
               PYNUCR = PYNUCR + PYNUCL (IKPMX)
               PZNUCR = PZNUCR + PZNUCL (IKPMX)
               ENUCR  = ENUCR  + ENNUC  (IKPMX)
               IBNUCR = IBNUCR + IBAR   (KPNUC)
               ICNUCR = ICNUCR + ICH    (KPNUC)
            END IF
            BNPREV = BNPREV + BNDUSE
            IF ( IREINT .LE. 0 ) THEN
               LPRCYC = .TRUE.
            ELSE
               LPRCYC = .FALSE.
            END IF
            GO TO 200
         ELSE IF ( BIMPCT .GT. RADTOT ) THEN
            KRFNUC (IKPMX) = KRFNUC (IKPMX) + 1
            CALL GRNDM(RNDM,1)
            SINTHE = RNDM ( 1 )
            COSTHE = SQRT ( 1.D+00 - SINTHE )
            SINTHE = SQRT ( SINTHE )
  400       CONTINUE
               CALL GRNDM(RNDM,2)
               RPHI1 = 2.D+00 * RNDM (1) - 1.D+00
               RPHI2 = 2.D+00 * RNDM (2) - 1.D+00
               RPHI12 = RPHI1 * RPHI1
               RPHI22 = RPHI2 * RPHI2
               RSQ = RPHI12 + RPHI22
            IF ( RSQ .GT. 1.D+00 ) GO TO 400
            SINPHI = 2.D+00 * RPHI1 * RPHI2 / RSQ
            COSPHI = ( RPHI12 - RPHI22 ) / RSQ
            SINT02 = CXIMPC**2 + CYIMPC**2
            IF ( SINT02 .LT. ANGLSQ ) THEN
               PXNUCL (IKPMX) = COSPHI * SINTHE * PNUCCO
               PYNUCL (IKPMX) = SINPHI * SINTHE * PNUCCO
               PZNUCL (IKPMX) = CZIMPC * COSTHE * PNUCCO
            ELSE
               SINTH0 = SQRT ( SINT02 )
               UPRIME = SINTHE * COSPHI
               VPRIME = SINTHE * SINPHI
               COSPH0 = CXIMPC / SINTH0
               SINPH0 = CYIMPC / SINTH0
               PXNUCL (IKPMX) = ( UPRIME * COSPH0 * CZIMPC - VPRIME
     &                        * SINPH0 + COSTHE * CXIMPC ) * PNUCCO
               PYNUCL (IKPMX) = ( UPRIME * SINPH0 * CZIMPC + VPRIME
     &                        * COSPH0 + COSTHE * CYIMPC ) * PNUCCO
               PZNUCL (IKPMX) = ( - UPRIME * SINTH0 + COSTHE * CZIMPC )
     &                        * PNUCCO
            END IF
            PNUCL  (IKPMX) = PNUCCO
            XSTNUC (IKPMX) = XIMPTR
            YSTNUC (IKPMX) = YIMPTR
            ZSTNUC (IKPMX) = ZIMPTR
            RSTNUC (IKPMX) = ABS (RIMPTR)
            ENNUC (IKPMX) = EKECON + AM (KPNUC)
            EKFNUC (IKPMX) = -AINFNT
            GO TO 200
         ELSE
            IF ( ( KPNUC .EQ. 1 .OR. KPNUC .EQ. 8 ) .AND. EKECON .LE.
     &           ETHREI ) THEN
               KPNUCL (IKPMX) = - KPNUCL (IKPMX)
               ENNUC (IKPMX) = EKECON + AM (KPNUC)
               NPROT = NPROT + ICH (KPNUC)
               NNEUT = NNEUT + 1 - ICH (KPNUC)
               IF ( .NOT. LBCHCK .AND. IKPNWI .GT. 0 ) THEN
                  LPRCYC = .TRUE.
                  NNUCTS = NNUCTS + 1
                  INUCTS (NNUCTS) = IKPMX
                  JNUCTS = NUSCIN
                  EKFNUC (IKPMX) = EKFPRE
                  ENNUC  (IKPMX) = ENNUC (IKPMX) - BNDUSE + BNDGEN
                  RSTNUC (IKPMX) = BNDGEN
                  LPRCYC = .TRUE.
               ELSE
                  EKFNUC (IKPMX) = -AINFNT
                  IREINT = IREINT + 1
                  IF ( EKFREI .LT. ANGLGB ) THEN
                     EKFREI = 0.5D+00 * ( EKFIMP + EKFPRO )
                     RHOREI = 0.5D+00 * ( RHOIMP + RHOIMT )
                  END IF
                  LPRCYC = .FALSE.
               END IF
               GO TO 200
            END IF
            LBIMPC = .FALSE.
            KPRIN  = KPNUC
            KPNUCL (IKPMX) = 0
            KRFLIN = KRFNUC (IKPMX)
            IF ( EKNNUC .NE. EKECON ) PNUCCO = SQRT ( EKNNUC * ( EKNNUC
     &         + 2.D+00 * AM (KPRIN) ) )
            CXIMPC = PXNUCL (IKPMX) / PNUCL (IKPMX)
            CYIMPC = PYNUCL (IKPMX) / PNUCL (IKPMX)
            CZIMPC = PZNUCL (IKPMX) / PNUCL (IKPMX)
            XSTNUC (IKPMX) = XIMPTR
            YSTNUC (IKPMX) = YIMPTR
            ZSTNUC (IKPMX) = ZIMPTR
            RSTNUC (IKPMX) = ABS (RIMPTR)
            GO TO 100
         END IF
      END IF
      NEXPEM = 0
      IF ( LGDHPR ) THEN
         LBCHCK = .TRUE.
         EKECON = EKE
         PNUCCO = PPROJ
         CALL BIMSEL ( KPROJ, TXX, TYY, TZZ, LBCHCK )
         LELSTC = .FALSE.
         RHOIMP = 0.5D+00 * ( RHOIMP + RHOIMT )
         EKFIMP = 0.5D+00 * ( EKFIMP + EKFPRO )
         RHOMEM = RHOIMP
         EKFMEM = EKFIMP
      END IF
*
      ANOW  = BBRES
      ZNOW  = ZZRES
      PXRES = PXORI
      PYRES = PYORI
      PZRES = PZORI
      PTRES = PTORI
      ERES  = EKE  + AM (KPROJ) + AMNTAR
      IF ( LGDHPR ) THEN
         IF ( KPROJ .EQ. 1 ) THEN
            NPROT = NPROT + 1
         ELSE IF ( KPROJ .EQ. 8 ) THEN
            NNEUT  = NNEUT + 1
         END IF
         ACOLL = BBTAR - 1.D+00
         IF ( KNUCIM .EQ. 1 ) THEN
            NPROT = NPROT + 1
            ZCOLL = ZZTAR - 1.D+00
         ELSE
            NNEUT = NNEUT + 1
            ZCOLL = ZZTAR
         END IF
         NHOLE = NHOLE + 1
      ELSE
         ACOLL = BBTAR - 1.D+00
         IF ( KPROJ .EQ. 1 ) THEN
            NPROT = NPROT + 1
            PRPONP = ZNOW / ( 3.D+00 * ANOW - 2.D+00 * ZNOW )
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. PRPONP ) THEN
               NPROT = NPROT + 1
               ZCOLL = ZZTAR - 1.D+00
               IPRTYP = KPROJ * 10 + 1
            ELSE
               NNEUT = NNEUT + 1
               ZCOLL = ZZTAR
               IPRTYP = KPROJ * 10 + 8
            END IF
            NHOLE = NHOLE + 1
         ELSE IF ( KPROJ .EQ. 8 ) THEN
            NNEUT  = NNEUT + 1
            PRNONP = 3.D+00 * ZNOW / ( 2.D+00 * ZNOW + ANOW )
            CALL GRNDM(RNDM,1)
            IF ( RNDM (1) .LT. PRNONP ) THEN
               NPROT = NPROT + 1
               ZCOLL = ZZTAR  - 1.D+00
               IPRTYP = KPROJ * 10 + 1
            ELSE
               NNEUT = NNEUT + 1
               ZCOLL = ZZTAR
               IPRTYP = KPROJ * 10 + 8
            END IF
            NHOLE = NHOLE + 1
         ELSE
            STOP 'KPROJ'
         END IF
      END IF
      ICYCL  = 0
      LEMISS = .FALSE.
  500 CONTINUE
      CALL PREPRE ( WEE, NNEUT, NPROT, NHOLE, ICYCL )
  530 CONTINUE
      IF ( IBRES .GT. 0 ) THEN
         EKR0   = ERES - AMNRES
         ATTNUM = ELBNDE (ICHTAR) - ELBNDE (ICRES) - EKR0 * ( AMMRES
     &          - AMNRES ) / AMMRES
         ERES   = ERES + AMMTAR - AMNTAR - ( ZZTAR - ZNOW ) * AMELEC
     &          + ATTNUM
         EKRES  = ERES - AMMRES
      ELSE
         AMMRES = 0.D+00
         AMNRES = 0.D+00
         ERES   = 0.D+00
         EKR0   = 0.D+00
         EKRES  = 0.D+00
         TVTENT = 0.D+00
         GO TO 600
      END IF
      IF ( EKRES .LE. 0.D+00 ) THEN
         WRITE ( LUNERR,* )' Peanut: negative kinetic energy for',
     &                     ' the residual nucleus!!',ICRES,IBRES,
     &                       REAL (EKRES)
         IF ( EKRES .LT. -3.D-3 ) THEN
            LRESMP = .TRUE.
            RETURN
         END IF
         EKRES  = 0.D+00
         TVRECL = 0.D+00
         AMSTAR = AMMRES
         TVCMS  = 0.D+00
         PTRES2 = 0.D+00
         PXRES  = 0.D+00
         PYRES  = 0.D+00
         PZRES  = 0.D+00
         PTRES  = 0.D+00
      ELSE
         PTRES2 = PTRES * PTRES
         AMSTAR = ( ERES - PTRES ) * ( ERES + PTRES )
         IF ( AMSTAR .GE. AMMRES**2 ) THEN
            AMSTAR = SQRT ( AMSTAR )
            TVCMS  = AMSTAR - AMMRES
         ELSE IF ( AMMRES**2 - AMSTAR .LT. 2.D+00 * AMSTAR * TVEPSI
     &             ) THEN
            AMSTAR = AMMRES
            ERES   = SQRT ( AMSTAR**2 + PTRES**2 )
            TVCMS  = 0.D+00
         ELSE IF ( AMSTAR .LE. 0.D+00 ) THEN
            WRITE ( LUNERR,* )' Peanut: immaginary mass for',
     &                        ' the residual nucleus!!',ICRES,IBRES,
     &                          REAL (AMSTAR)
            LRESMP = .TRUE.
            RETURN
         ELSE
            AMSTAR = SQRT ( AMSTAR )
            IF ( AMMRES - AMSTAR .LT. TVEPSI ) THEN
               AMSTAR = AMMRES
               TVCMS  = 0.D+00
               TVRECL = ERES - AMSTAR
               GO TO 550
            END IF
            WRITE ( LUNERR,* )' Peanut: negative excitation energy for',
     &                        ' the residual nucleus!!',ICRES,IBRES,
     &                          REAL ( AMSTAR - AMMRES )
            IF ( AMSTAR - AMMRES .LT. -3.D-3 ) THEN
               LRESMP = .TRUE.
               RETURN
            END IF
            AMSTAR = AMMRES
            TVCMS  = 0.D+00
            HELP   = SQRT ( ( ERES - AMMRES ) * ( ERES + AMMRES ) )
     &             / PTRES
            PXRES = PXRES * HELP
            PYRES = PYRES * HELP
            PZRES = PZRES * HELP
            PTRES = PTRES * HELP
         END IF
         TVRECL = ERES - AMSTAR
      END IF
  550 CONTINUE
      IF ( TVRECL .LT. 0.D+00 ) THEN
         TVRECL = 0.D+00
      END IF
      TV     = 0.D+00
      EKRES  = TVRECL
  600 CONTINUE
      EOTEST = EOTEST - EINTR - ENUCR - EKR0 - AMNRES
      IF ( ABS (EOTEST) .GT. ETEPS ) THEN
         WRITE (LUNERR,*)' Peanut: eotest failure',EOTEST
         LRESMP = .TRUE.
         RETURN
      END IF
      IF ( IBRES .EQ. 0 ) RETURN
      EOTEST = ETTOT + AMMTAR - AMNTAR + ATTNUM
      IF ( KPROJ .EQ. 1 ) THEN
         EOTEST = EOTEST + AMHEAV (2) - AM (1)
      ELSE IF ( KPROJ .EQ. 8 ) THEN
         EOTEST = EOTEST + AMHEAV (1) - AM (8)
      ELSE
         EOTEST = EOTEST + ICH(KPROJ) * AMELEC
      END IF
      IF ( LEVPRT ) THEN
         CALL EVEVAP ( WEE )
         IF ( LRESMP ) RETURN
      ELSE
         TVHEAV = 0.D+00
      END IF
      DO 1000 KP = NP0+1,NP
         IF ( KPART (KP) .EQ. 1 ) THEN
            EOTEST = EOTEST - TKI (KP) - AMHEAV (2)
            IBTOT  = IBTOT  - 1
            ICHTOT = ICHTOT - 1
         ELSE IF ( KPART (KP) .EQ. 8 ) THEN
            EOTEST = EOTEST - TKI (KP) - AMHEAV (1)
            IBTOT  = IBTOT  - 1
         ELSE
            EOTEST = EOTEST - TKI (KP) - AM (KPART(KP))
            IBTOT  = IBTOT  - IBAR (KPART(KP))
            ICHTOT = ICHTOT - ICH  (KPART(KP))
         END IF
 1000 CONTINUE
      EOTEST = EOTEST - TVHEAV - IEVDEU * AMHEAV (3)
     &       - IEVTRI * AMHEAV (4)
     &       - IEV3HE * AMHEAV (5)
     &       - IEV4HE * AMHEAV (6)
     &       - AMMRES - TVRECL
      IBTOT  = IBTOT  - IEVDEU * 2 - IEVTRI * 3 - IEV3HE * 3
     &       - IEV4HE * 4
      ICHTOT = ICHTOT - IEVDEU - IEVTRI - IEV3HE * 2
     &       - IEV4HE * 2
      IF ( LRNFSS ) THEN
         IF ( LHEAVY ) THEN
            DO 2000 JP = 1, NPHEAV
               IF ( KHEAVY (JP) .GT. 6 ) THEN
                  EOTEST = EOTEST - AMHEAV (JP)
                  IBTOT  = IBTOT  - IBHEAV (KHEAVY(JP))
                  ICHTOT = ICHTOT - ICHEAV (KHEAVY(JP))
               END IF
 2000       CONTINUE
         ELSE
            DO 2100 JFISS = 1, NFISS
               IBHLP = NINT (ATFIS(JFISS))
               IF ( IBHLP .GT. 0 ) THEN
                  ICHLP  = NINT (ZTFIS(JFISS))
                  EOTEST = EOTEST - 1.D-03 * AMTFIS (JFISS)
                  IBTOT  = IBTOT  - IBHLP
                  ICHTOT = ICHTOT - ICHLP
               END IF
 2100       CONTINUE
         END IF
      END IF
      IF ( ABS (EOTEST) .GT. 1.D+3 * ETEPS ) THEN
         WRITE (LUNERR,*)
     &   ' Peanut failure!!, Eotest,Ammres,Tvrecl,Ibres,Icres',
     &                       EOTEST,AMMRES,TVRECL,IBRES,ICRES
      END IF
      IF ( IBTOT .NE. IBRES .OR. ICHTOT .NE. ICRES ) THEN
         WRITE (LUNERR,*)
     &   ' Peanut failure!!, Ichtot, Icres, Ibtot, Ibres',
     &                       ICHTOT, ICRES, IBTOT, IBRES
      END IF
*=== End of subroutine peanut =========================================*
      RETURN
      END
