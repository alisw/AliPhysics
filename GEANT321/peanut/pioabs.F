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
*$ CREATE PIOABS.FOR
*COPY PIOABS
*
*=== pioabs ===========================================================*
*
      SUBROUTINE PIOABS ( IKPMX , KRFLIN, WEE   , ERECMN, LBIMPC,
     &                    LBCHCK, ICYCL , NHOLE , NPROT , NNEUT ,
     &                    LEXIT , LNWINT )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*
#include "geant321/balanc.inc"
#include "geant321/finuc.inc"
#include "geant321/nucdat.inc"
#include "geant321/nucgeo.inc"
#include "geant321/parevt.inc"
#include "geant321/parnuc.inc"
#include "geant321/part.inc"
#include "geant321/resnuc.inc"
*
      COMMON / FKPLOC / IABCOU
      REAL RNDM(1)
      LOGICAL LBCHCK, LBIMPC, LTROUB, LEXIT, LNWINT
*
      NPNCLD = NPNUC
 1000 CONTINUE
      IF ( LABRST .OR. LABSRP ) THEN
         LNWINT = .FALSE.
      ELSE
         LEXIT  = .FALSE.
         LNWINT = .TRUE.
         RETURN
      END IF
      NHOLE  = NHOLE + NTARGT
      NTARLD = NTARGT
      ICYCL  = ICYCL + 1
      IF ( NTARGT .EQ. 1 ) THEN
         IF ( .NOT. LABRST )  STOP '???_rad_flight_abs'
         LABRST = .FALSE.
         LABSRP = .FALSE.
         PFROUT = PFRIMP / PFRCEN (1) * PFRCEN (2)
         EKFOUT = SQRT ( AMNUSQ (2) + PFROUT**2 ) - AMNUCL (2)
         POTINC = EKEWLL - EKECON + EKFERM
         POTOUT = EKFERM + EKFOUT + BNENRG (2) - EKFIMP - BNENRG (1)
         ERES   = EKEWLL + AM (KPRIN) + EKFERM + AM (KNUCIM)
     &          + POTOUT - POTINC
         AMNREC = AMNTAR - AMUC12
         ERECMN = 0.5D+00 * PTRES2 / AMNREC
         ERECMN = 0.D+00
         UMO2   = ERES*ERES - PTRES2
         UMO    = SQRT (UMO2)
         GAMCM = ERES  / UMO
         ETAX  = PXRES / UMO
         ETAY  = PYRES / UMO
         ETAZ  = PZRES / UMO
         ECMSNU = 0.5D+00 * ( UMO2 + AMNUSQ (2) ) / UMO
         PCMS   = UMO - ECMSNU
         CALL RACO ( PCMSX, PCMSY, PCMSZ )
         PCMSX = PCMS * PCMSX
         PCMSY = PCMS * PCMSY
         PCMSZ = PCMS * PCMSZ
         NPNUC = NPNUC + 1
         KPNUCL (NPNUC) = 7
         KRFNUC (NPNUC) = KRFLIN + 1
         ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
         PHELP  = PCMS + ETAPCM / ( GAMCM + 1.D+00 )
         ENNUC  (NPNUC) = GAMCM * PCMS + ETAPCM
         PXHELP = PCMSX + ETAX * PHELP
         PYHELP = PCMSY + ETAY * PHELP
         PZHELP = PCMSZ + ETAZ * PHELP
         PXRES = PXRES - PXHELP
         PYRES = PYRES - PYHELP
         PZRES = PZRES - PZHELP
         ERES  = ERES  - ENNUC (NPNUC)
         PTRES2= PXRES**2 + PYRES**2 + PZRES**2
         PXHLP = PXTTOT - PXHELP
         PYHLP = PYTTOT - PYHELP
         PZHLP = PZTTOT - PZHELP
         UMO2  = ( ETTOT - ENNUC (NPNUC) )**2 - PXHLP**2 - PYHLP**2
     &         - PZHLP**2
         EEXMNM = 0.D+00
         DELTU2 = UMO2 - ( AMNRES + EEXMNM )**2
         IF ( DELTU2 .LT. 0.D+00 ) THEN
            NPNUC  = NPNUC - 1
            LBCHCK = .FALSE.
            IF ( LBIMPC ) THEN
               CALL BIMNXT ( LBCHCK )
               RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
               EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
            ELSE
               CALL NWINXT ( LBCHCK )
               IF ( BIMPCT .GT. RADTOT ) THEN
                  NHOLE = NHOLE - NTARLD
                  ICYCL = ICYCL - 1
                  CALL PHDSET ( IKPMX )
                  IBRES = IBRES - IBAR (KPRIN)
                  ICRES = ICRES - ICH  (KPRIN)
                  BBRES = IBRES
                  ZZRES = ICRES
                  AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                   ( BBRES, ZZRES)
                  AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
                  LTROUB = .FALSE.
                  CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
                  IF ( LTROUB ) THEN
                     KPNUCL (IKPMX) = 0
                     UMO2  = ERES**2 - PTRES2
                     UMO = SQRT (UMO2)
                     WRITE ( LUNOUT,* )' PIO0_P:UMO,AMNRES',UMO,AMNRES
                     LEXIT = .TRUE.
                     RETURN
                  END IF
                  NP = NP + 1
                  TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
                  KPART (NP) = KPRIN
                  PLR   (NP) = PNUCL  (IKPMX)
                  CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
                  CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
                  CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
                  WEI   (NP) = WEE
                  KPNUCL (IKPMX) = 0
                  IOTHER = IOTHER + 1
                  PXNUCR = PXNUCR + PXNUCL (IKPMX)
                  PYNUCR = PYNUCR + PYNUCL (IKPMX)
                  PZNUCR = PZNUCR + PZNUCL (IKPMX)
                  ENUCR  = ENUCR  + ENNUC  (IKPMX)
                  IBNUCR = IBNUCR + IBAR   (KPART(NP))
                  ICNUCR = ICNUCR + ICH    (KPART(NP))
                  LEXIT  = .TRUE.
                  RETURN
               END IF
               XSTNUC (IKPMX) = XIMPTR
               YSTNUC (IKPMX) = YIMPTR
               ZSTNUC (IKPMX) = ZIMPTR
               RSTNUC (IKPMX) = ABS (RIMPTR)
            END IF
            NHOLE = NHOLE - NTARLD
            ICYCL = ICYCL - 1
            GO TO 1000
         END IF
         EKFNUC (NPNUC) = - AINFNT
         PXNUCL (NPNUC) = PXHELP
         PYNUCL (NPNUC) = PYHELP
         PZNUCL (NPNUC) = PZHELP
         PNUCL  (NPNUC) = ENNUC (NPNUC)
         XSTNUC (NPNUC) = XIMPTR
         YSTNUC (NPNUC) = YIMPTR
         ZSTNUC (NPNUC) = ZIMPTR
         RSTNUC (NPNUC) = ABS (RIMPTR)
         NPNUC = NPNUC + 1
         KPNUCL (NPNUC) = 8
         KRFNUC (NPNUC) = KRFLIN + 1
         ETAPCM = - ETAPCM
         PHELP  = ECMSNU + ETAPCM / ( GAMCM + 1.D+00 )
         ENNUC  (NPNUC) = GAMCM * ECMSNU + ETAPCM
         DEFRNU = DEFNEU
         IF ( ENNUC (NPNUC) - AM (8) .LE. EKFOUT + DEFRNU ) THEN
            NPNUC  = NPNUC - 2
            LBCHCK = .FALSE.
            IF ( LBIMPC ) THEN
               CALL BIMNXT ( LBCHCK )
               RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
               EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
            ELSE
               CALL NWINXT ( LBCHCK )
               IF ( BIMPCT .GT. RADTOT ) THEN
                  NHOLE = NHOLE - NTARLD
                  ICYCL = ICYCL - 1
                  CALL PHDSET ( IKPMX )
                  IBRES = IBRES - IBAR (KPRIN)
                  ICRES = ICRES - ICH  (KPRIN)
                  BBRES = IBRES
                  ZZRES = ICRES
                  AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                   ( BBRES, ZZRES)
                  AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
                  LTROUB = .FALSE.
                  CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
                  IF ( LTROUB ) THEN
                     KPNUCL (IKPMX) = 0
                     UMO2  = ERES**2 - PTRES2
                     UMO = SQRT (UMO2)
                     WRITE ( LUNOUT,* )' PIO0_T:UMO,AMNRES',UMO,AMNRES
                     LEXIT = .TRUE.
                     RETURN
                  END IF
                  NP = NP + 1
                  TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
                  KPART (NP) = KPRIN
                  PLR   (NP) = PNUCL  (IKPMX)
                  CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
                  CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
                  CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
                  WEI   (NP) = WEE
                  KPNUCL (IKPMX) = 0
                  IOTHER = IOTHER + 1
                  PXNUCR = PXNUCR + PXNUCL (IKPMX)
                  PYNUCR = PYNUCR + PYNUCL (IKPMX)
                  PZNUCR = PZNUCR + PZNUCL (IKPMX)
                  ENUCR  = ENUCR  + ENNUC  (IKPMX)
                  IBNUCR = IBNUCR + IBAR   (KPART(NP))
                  ICNUCR = ICNUCR + ICH    (KPART(NP))
                  LEXIT  = .TRUE.
                  RETURN
               END IF
               XSTNUC (IKPMX) = XIMPTR
               YSTNUC (IKPMX) = YIMPTR
               ZSTNUC (IKPMX) = ZIMPTR
               RSTNUC (IKPMX) = ABS (RIMPTR)
            END IF
            NHOLE = NHOLE - NTARLD
            ICYCL = ICYCL - 1
            GO TO 1000
         END IF
         EKFNUC (NPNUC) = EKFOUT
         PXNUCL (NPNUC) = -PCMSX + ETAX * PHELP
         PYNUCL (NPNUC) = -PCMSY + ETAY * PHELP
         PZNUCL (NPNUC) = -PCMSZ + ETAZ * PHELP
         PNUCL  (NPNUC) = SQRT ( PXNUCL (NPNUC)**2 + PYNUCL (NPNUC)**2
     &                         + PZNUCL (NPNUC)**2 )
         XSTNUC (NPNUC) = XIMPCT
         YSTNUC (NPNUC) = YIMPCT
         ZSTNUC (NPNUC) = ZIMPCT
         RSTNUC (NPNUC) = ABS (RIMPCT)
         LBIMPC = .FALSE.
         LEXIT  = .FALSE.
         NUSCIN = NUSCIN + 1
         ISCTYP (NUSCIN) = - ( KPRIN * 100 + KNUCIM )
         IF ( NUSCIN .EQ. 1 ) IPRTYP = ISCTYP (1)
         NHLEXP = NHLEXP + 1
         HOLEXP (NHLEXP) = EKFIMP - EKFERM
         RHOACT = 0.5D+00 * ( RHOIMP + RHOIMT )
         RHOEXP = RHOEXP + RHOACT
         EKFEXP = EKFEXP + 0.5D+00 * ( EKFIMP + EKFPRO )
         CALL NCLVFX
      ELSE
         LABRST = .FALSE.
         LABSRP = .FALSE.
         ITFRMI = 1 + KNUCIM / 8
         ITFRM2 = 1 + KNUCI2 / 8
         IF ( ICH (KPRIN) .GT. 0 ) THEN
            IOFRMI = 1
            IOFRM2 = ITFRM2
            DEFRPR = DEFPRO
            IF ( IOFRMI .EQ. ITFRM2 ) THEN
               EKFOUT = EKFIM2
               DEFRNU = DEFPRO
            ELSE
               PFROUT = PFRIMP / PFRCEN (1) * PFRCEN (2)
               EKFOUT = SQRT ( AMNUSQ (1) + PFROUT**2 ) - AMNUCL (1)
               DEFRNU = DEFNEU
            END IF
         ELSE IF ( ICH (KPRIN) .LT. 0 ) THEN
            IOFRMI = 2
            IOFRM2 = ITFRM2
            DEFRPR = DEFNEU
            IF ( IOFRMI .EQ. ITFRM2 ) THEN
               EKFOUT = EKFIM2
               DEFRNU = DEFNEU
            ELSE
               PFROUT = PFRIMP / PFRCEN (2) * PFRCEN (1)
               EKFOUT = SQRT ( AMNUSQ (2) + PFROUT**2 ) - AMNUCL (2)
               DEFRNU = DEFPRO
            END IF
         ELSE
            IOFRMI = ITFRMI
            IOFRM2 = ITFRM2
            EKFOUT = EKFIMP
            IF ( ITFRMI .EQ. 1 ) THEN
               DEFRPR = DEFPRO
            ELSE
               DEFRPR = DEFNEU
            END IF
            IF ( ITFRM2 .EQ. 1 ) THEN
               DEFRNU = DEFPRO
            ELSE
               DEFRNU = DEFNEU
            END IF
         END IF
         POTINC = EKEWLL - EKECON + EKFERM + EKFER2
         POTOUT = EKFERM + EKFER2 + EKFOUT + BNENRG (IOFRMI) - EKFIMP
     &          - BNENRG (ITFRMI)
         ERES   = EKEWLL + AM (KPRIN) + EKFERM + AM (KNUCIM)
     &          + EKFER2 + AM (KNUCI2) + POTOUT - POTINC
         AMNREC = AMNTAR - 2.D+00 * AMUC12
         PHLPSQ = ( PXRES - CXIMPC * PNUCCO )**2
     &          + ( PYRES - CYIMPC * PNUCCO )**2
     &          + ( PZRES - CZIMPC * PNUCCO )**2
         ERECMN = 0.5D+00 * PHLPSQ / AMNREC**2
         ERECMN = AMNREC * ERECMN * ( 1.D+00 - 0.25D+00 * ERECMN )
         ERECMN = 0.D+00
         UMO2   = ERES*ERES - PTRES2
         UMO    = SQRT (UMO2)
         GAMCM = ERES  / UMO
         ETAX  = PXRES / UMO
         ETAY  = PYRES / UMO
         ETAZ  = PZRES / UMO
         ECMSPR = 0.5D+00 * ( UMO2 + AMNUSQ (IOFRMI) - AMNUSQ (IOFRM2) )
     &          / UMO
         ECMSNU = UMO - ECMSPR
         PCMS   = SQRT ( ( ECMSPR - AMNUCL (IOFRMI) ) * ( ECMSPR
     &          + AMNUCL (IOFRMI) ) )
         CALL RACO ( PCMSX, PCMSY, PCMSZ )
         PCMSX = PCMS * PCMSX
         PCMSY = PCMS * PCMSY
         PCMSZ = PCMS * PCMSZ
         NPNUC = NPNUC + 1
         KPNUCL (NPNUC) = 1 + 7 * ( IOFRMI - 1 )
         KRFNUC (NPNUC) = KRFLIN + 1
         ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
         PHELP  = ECMSPR + ETAPCM / ( GAMCM + 1.D+00 )
         ENNUC  (NPNUC) = GAMCM * ECMSPR + ETAPCM
         IF ( ENNUC (NPNUC) - AMNUCL (IOFRMI) .LE. EKFOUT + DEFRPR )THEN
            NPNUC  = NPNUC - 1
            LBCHCK = .FALSE.
            IF ( LBIMPC ) THEN
               CALL BIMNXT ( LBCHCK )
               RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
               EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
            ELSE
               CALL NWINXT ( LBCHCK )
               IF ( BIMPCT .GT. RADTOT ) THEN
                  NHOLE = NHOLE - NTARLD
                  ICYCL = ICYCL - 1
                  CALL PHDSET ( IKPMX )
                  IBRES = IBRES - IBAR (KPRIN)
                  ICRES = ICRES - ICH  (KPRIN)
                  BBRES = IBRES
                  ZZRES = ICRES
                  AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                   ( BBRES, ZZRES)
                  AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
                  LTROUB = .FALSE.
                  CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
                  IF ( LTROUB ) THEN
                     KPNUCL (IKPMX) = 0
                     UMO2  = ERES**2 - PTRES2
                     UMO = SQRT (UMO2)
                     WRITE ( LUNOUT,* )' PIO0_P:UMO,AMNRES',UMO,AMNRES
                     LEXIT = .TRUE.
                     RETURN
                  END IF
                  NP = NP + 1
                  TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
                  KPART (NP) = KPRIN
                  PLR   (NP) = PNUCL  (IKPMX)
                  CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
                  CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
                  CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
                  WEI   (NP) = WEE
                  KPNUCL (IKPMX) = 0
                  IOTHER = IOTHER + 1
                  PXNUCR = PXNUCR + PXNUCL (IKPMX)
                  PYNUCR = PYNUCR + PYNUCL (IKPMX)
                  PZNUCR = PZNUCR + PZNUCL (IKPMX)
                  ENUCR  = ENUCR  + ENNUC  (IKPMX)
                  IBNUCR = IBNUCR + IBAR   (KPART(NP))
                  ICNUCR = ICNUCR + ICH    (KPART(NP))
                  LEXIT  = .TRUE.
                  RETURN
               END IF
               XSTNUC (IKPMX) = XIMPTR
               YSTNUC (IKPMX) = YIMPTR
               ZSTNUC (IKPMX) = ZIMPTR
               RSTNUC (IKPMX) = ABS (RIMPTR)
            END IF
            NHOLE = NHOLE - NTARLD
            ICYCL = ICYCL - 1
            GO TO 1000
         END IF
         EKFNUC (NPNUC) = EKFOUT
         PXNUCL (NPNUC) = PCMSX + ETAX * PHELP
         PYNUCL (NPNUC) = PCMSY + ETAY * PHELP
         PZNUCL (NPNUC) = PCMSZ + ETAZ * PHELP
         PNUCL  (NPNUC) = SQRT ( PXNUCL (NPNUC)**2 + PYNUCL (NPNUC)**2
     &                         + PZNUCL (NPNUC)**2 )
         XSTNUC (NPNUC) = XIMPTR
         YSTNUC (NPNUC) = YIMPTR
         ZSTNUC (NPNUC) = ZIMPTR
         RSTNUC (NPNUC) = ABS (RIMPTR)
         NPNUC = NPNUC + 1
         KPNUCL (NPNUC) = 1 + 7 * ( IOFRM2 - 1 )
         KRFNUC (NPNUC) = KRFLIN + 1
         ETAPCM = - ETAPCM
         PHELP  = ECMSNU + ETAPCM / ( GAMCM + 1.D+00 )
         ENNUC  (NPNUC) = GAMCM * ECMSNU + ETAPCM
         IF ( ENNUC (NPNUC) - AMNUCL (IOFRM2) .LE. EKFIM2 + DEFRNU )THEN
            NPNUC  = NPNUC - 2
            LBCHCK = .FALSE.
            IF ( LBIMPC ) THEN
               CALL BIMNXT ( LBCHCK )
               RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
               EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
            ELSE
               CALL NWINXT ( LBCHCK )
               IF ( BIMPCT .GT. RADTOT ) THEN
                  NHOLE = NHOLE - NTARLD
                  ICYCL = ICYCL - 1
                  CALL PHDSET ( IKPMX )
                  IBRES = IBRES - IBAR (KPRIN)
                  ICRES = ICRES - ICH  (KPRIN)
                  BBRES = IBRES
                  ZZRES = ICRES
                  AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                   ( BBRES, ZZRES)
                  AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
                  LTROUB = .FALSE.
                  CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
                  IF ( LTROUB ) THEN
                     KPNUCL (IKPMX) = 0
                     UMO2  = ERES**2 - PTRES2
                     UMO = SQRT (UMO2)
                     WRITE ( LUNOUT,* )' PIO0_T:UMO,AMNRES',UMO,AMNRES
                     LEXIT = .TRUE.
                     RETURN
                  END IF
                  NP = NP + 1
                  TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
                  KPART (NP) = KPRIN
                  PLR   (NP) = PNUCL  (IKPMX)
                  CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
                  CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
                  CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
                  WEI   (NP) = WEE
                  KPNUCL (IKPMX) = 0
                  IOTHER = IOTHER + 1
                  PXNUCR = PXNUCR + PXNUCL (IKPMX)
                  PYNUCR = PYNUCR + PYNUCL (IKPMX)
                  PZNUCR = PZNUCR + PZNUCL (IKPMX)
                  ENUCR  = ENUCR  + ENNUC  (IKPMX)
                  IBNUCR = IBNUCR + IBAR   (KPART(NP))
                  ICNUCR = ICNUCR + ICH    (KPART(NP))
                  LEXIT  = .TRUE.
                  RETURN
               END IF
               XSTNUC (IKPMX) = XIMPTR
               YSTNUC (IKPMX) = YIMPTR
               ZSTNUC (IKPMX) = ZIMPTR
               RSTNUC (IKPMX) = ABS (RIMPTR)
            END IF
            NHOLE = NHOLE - NTARLD
            ICYCL = ICYCL - 1
            GO TO 1000
         END IF
         EKFNUC (NPNUC) = EKFIM2
         PXNUCL (NPNUC) = -PCMSX + ETAX * PHELP
         PYNUCL (NPNUC) = -PCMSY + ETAY * PHELP
         PZNUCL (NPNUC) = -PCMSZ + ETAZ * PHELP
         PNUCL  (NPNUC) = SQRT ( PXNUCL (NPNUC)**2 + PYNUCL (NPNUC)**2
     &                         + PZNUCL (NPNUC)**2 )
         XSTNUC (NPNUC) = XIMPCT
         YSTNUC (NPNUC) = YIMPCT
         ZSTNUC (NPNUC) = ZIMPCT
         RSTNUC (NPNUC) = ABS (RIMPCT)
         LBIMPC = .FALSE.
         LEXIT  = .FALSE.
         NUSCIN = NUSCIN + 1
         ISCTYP (NUSCIN) = - ( KPRIN * 100 + KNUCIM * 10 + KNUCI2 )
         IF ( NUSCIN .EQ. 1 ) IPRTYP = ISCTYP (1)
         IABCOU = IABCOU + 1
         NHLEXP = NHLEXP + 2
         HOLEXP (NHLEXP-1) = EKFIMP - EKFERM
         HOLEXP (NHLEXP)   = EKFIM2 - EKFER2
         RHOACT = 0.6666666666666666D+00 * RHOIMP
     &          + 0.3333333333333333D+00 * RHOIMT
         RHOEXP = RHOEXP + 2.D+00 * RHOACT
         EKFEXP = EKFEXP + 0.6666666666666666D+00 * ( EKFIMP + EKFIM2
     &          + EKFPRO )
         CALL NCLVFX
      END IF
      DO 3000 KP = NPNCLD+1, NPNUC
         KPNUC = KPNUCL (KP)
         IF ( AM (KPNUC) .LE. 0.D+00 ) THEN
            TAUTAU = RZNUCL / PNUCL (KP)
         ELSE
            TAUEFF = 0.5D+00 * TAUFOR * AM (13) / AM (KPNUC)
            CALL GRNDM(RNDM,1)
            TAUTAU = - TAUEFF / AM (KPNUC) * LOG ( 1.D+00 - RNDM
     &             (1) )
            TAUTAU = MAX ( TAUTAU, RZNUCL / PNUCL (KP) )
         END IF
         XSTNUC (KP) = XSTNUC (KP) + PXNUCL (KP) * TAUTAU
         YSTNUC (KP) = YSTNUC (KP) + PYNUCL (KP) * TAUTAU
         ZSTNUC (KP) = ZSTNUC (KP) + PZNUCL (KP) * TAUTAU
         RSTNUC (KP) = SQRT ( XSTNUC (KP)**2 + YSTNUC (KP)**2
     &               + ZSTNUC (KP)**2 )
         RHNUCL (KP) = RHOACT
 3000 CONTINUE
      RETURN
*=== End of subroutine pioabs =========================================*
      END
