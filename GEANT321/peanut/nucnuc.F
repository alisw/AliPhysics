*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:22:01  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.46  by  S.Giani
*-- Author :
*$ CREATE NUCNUC.FOR
*COPY NUCNUC
*
*=== nucnuc ===========================================================*
*
      SUBROUTINE NUCNUC ( IKPMX , KRFLIN, WEE   , ERECMN, LBIMPC,
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
      REAL RNDM(2)
      LOGICAL LBCHCK, LBIMPC, LTROUB, LEXIT, LNWINT
*
      NPNCLD = NPNUC
 1000 CONTINUE
      IF ( LELSTC ) THEN
         LELSTC = .FALSE.
         LNWINT = .FALSE.
      ELSE
         LEXIT  = .FALSE.
         LNWINT = .TRUE.
         RETURN
      END IF
      NHOLE = NHOLE + 1
      ICYCL = ICYCL + 1
      IF ( NUSCIN .LE. 0 ) THEN
         PXHLP  = PXFERM + PXPROJ - CXIMPC * PPRWLL
         PYHLP  = PYFERM + PYPROJ - CYIMPC * PPRWLL
         PZHLP  = PZFERM + PZPROJ - CZIMPC * PPRWLL
         ERECMN = 0.5D+00 * ( PXHLP**2 + PYHLP**2 + PZHLP**2 ) / AMNTAR
         ERECMN = ERECMN * ( 1.D+00 - 0.25D+00 * ERECMN )
         ERECMN = MAX ( ERECMN, EKFERM * AM (KNUCIM) / AMNTAR )
      ELSE
         ERECMN = EKFERM * AM (KNUCIM) / AMNTAR
      END IF
      ERECMN = 0.D+00
      POTINC = EKEWLL - EKECON + EKFERM
      IF ( NUSCIN .EQ. 0 .AND. KPRIN .NE. KNUCIM ) THEN
         ZAFT   = ZNOW + ICH (KPRIN) - ICH (KNUCIM)
         AMMAFT = ANOW * AMUAMU + 0.001D+00 * FKENER ( ANOW, ZAFT )
         AMNAFT = AMMAFT - ZAFT * AMELEC + ELBNDE (NINT(ZAFT))
         IF ( KPRIN .EQ. 1 ) THEN
            DEFRPR = DEFPRO
         ELSE
            DEFRPR = DEFNEU
         END IF
         IF ( KNUCIM .EQ. 8 ) THEN
            DEFRNU = MAX ( DEFNEU, HLFHLF * EEXANY + ERECMN
     &             + EKFERM - EKFIMP )
         ELSE
            DEFRNU = MAX ( DEFPRO, HLFHLF * EEXANY + ERECMN
     &             + EKFERM - EKFIMP )
         END IF
      ELSE IF ( NUSCIN .EQ. 0 ) THEN
         IF ( KPRIN .EQ. 1 ) THEN
            DEFRPR = MAX ( DEFPRO, HLFHLF * EEXANY + ERECMN
     &             + EKFERM - EKFIMP )
            DEFRNU = DEFRPR
         ELSE
            DEFRPR = MAX ( DEFNEU, HLFHLF * EEXANY + ERECMN
     &             + EKFERM - EKFIMP )
            DEFRNU = DEFRPR
         END IF
      ELSE
         ADDPRO = 0.D+00
         ADDTAR = 0.D+00
         IF ( KPRIN .EQ. 8 ) THEN
            DEFRPR = DEFNEU + ADDPRO
         ELSE
            DEFRPR = DEFPRO + ADDPRO
         END IF
         IF ( KNUCIM .EQ. 8 ) THEN
            DEFRNU = DEFNEU + ADDTAR
         ELSE
            DEFRNU = DEFPRO + ADDTAR
         END IF
      END IF
      NPLSIN = 2
      UMO2   = ERES**2 - PTRES2
      UMO    = SQRT (UMO2)
      GAMCM = ERES  / UMO
      ETAX  = PXRES / UMO
      ETAY  = PYRES / UMO
      ETAZ  = PZRES / UMO
      ETAPCM = ETAX * PXPROJ + ETAY * PYPROJ + ETAZ * PZPROJ
      ECMSPR = GAMCM * ( EKEWLL + AM (KPRIN) ) - ETAPCM
      PHELP  = EKEWLL + AM (KPRIN) - ETAPCM / ( GAMCM + 1.D+00 )
      PCMSX  = PXPROJ - ETAX * PHELP
      PCMSY  = PYPROJ - ETAY * PHELP
      PCMSZ  = PZPROJ - ETAZ * PHELP
      ETAPCM = ETAX * PXFERM + ETAY * PYFERM + ETAZ * PZFERM
      ECMSNU = GAMCM * ( EKFERM + AM (KNUCIM) ) - ETAPCM
      PCMS2  = PCMSX**2 + PCMSY**2 + PCMSZ**2
      PCMS   = SQRT ( PCMS2 )
      IF ( KPRIN .EQ. KNUCIM ) THEN
         EKESAM = 0.5D+00 * ( UMO2 - AM (KPRIN)**2 - AM (KNUCIM)**2 )
     &          / AM (KNUCIM) - AM (KPRIN)
         CALL SAMCST ( 1, EKESAM, COSTHE )
      ELSE IF ( KPRIN .EQ. 8 ) THEN
         EKESAM = 0.5D+00 * ( UMO2 - AM (KPRIN)**2 - AM (KNUCIM)**2 )
     &          / AM (KNUCIM) - AM (KPRIN)
         CALL SAMCST ( 8, EKESAM, COSTHE )
      ELSE
         EKESAM = 0.5D+00 * ( UMO2 - AM (KNUCIM)**2 - AM (KPRIN)**2 )
     &          / AM (KPRIN) - AM (KNUCIM)
         CALL SAMCST ( 8, EKESAM, COSTHE )
         COSTHE = - COSTHE
      END IF
      SINTHE = SQRT ( ( 1.D+00 - COSTHE ) * ( 1.D+00 + COSTHE ) )
 2000 CONTINUE
         CALL GRNDM(RNDM,2)
         RPHI1 = 2.D+00 * RNDM (1) - 1.D+00
         RPHI2 = 2.D+00 * RNDM (2) - 1.D+00
         RPHI12 = RPHI1 * RPHI1
         RPHI22 = RPHI2 * RPHI2
         RSQ = RPHI12 + RPHI22
      IF ( RSQ .GT. 1.D+00 ) GO TO 2000
      SINPHI = 2.D+00 * RPHI1 * RPHI2 / RSQ
      COSPHI = ( RPHI12 - RPHI22 ) / RSQ
      CZAXCM = PCMSZ / PCMS
      SINT02 = ( 1.D+00 - CZAXCM ) * ( 1.D+00 + CZAXCM )
      IF ( SINT02 .LT. ANGLSQ ) THEN
         CXCMS = COSPHI * SINTHE
         CYCMS = SINPHI * SINTHE
         CZCMS = CZAXCM * COSTHE
      ELSE
         SINTH0 = SQRT ( SINT02 )
         UPRIME = SINTHE * COSPHI
         VPRIME = SINTHE * SINPHI
         CXAXCM = PCMSX / PCMS
         CYAXCM = PCMSY / PCMS
         COSPH0 = CXAXCM / SINTH0
         SINPH0 = CYAXCM / SINTH0
         CXCMS = UPRIME * COSPH0 * CZAXCM - VPRIME * SINPH0
     &         + COSTHE * CXAXCM
         CYCMS = UPRIME * SINPH0 * CZAXCM + VPRIME * COSPH0
     &         + COSTHE * CYAXCM
         CZCMS = - UPRIME * SINTH0 + COSTHE * CZAXCM
      END IF
      PCMSX  = PCMS * CXCMS
      PCMSY  = PCMS * CYCMS
      PCMSZ  = PCMS * CZCMS
      NPNUC = NPNUC + 1
      KPNUCL (NPNUC) = KPRIN
      KRFNUC (NPNUC) = KRFLIN + 1
      ETAPCM = ETAX * PCMSX + ETAY * PCMSY + ETAZ * PCMSZ
      PHELP  = ECMSPR + ETAPCM / ( GAMCM + 1.D+00 )
      ENNUC  (NPNUC) = GAMCM * ECMSPR + ETAPCM
      IF ( ENNUC (NPNUC) - AM (KPRIN) .LE. EKFPRO + DEFRPR ) THEN
         NPNUC  = NPNUC - 1
         LBCHCK = .FALSE.
         IF ( LBIMPC ) THEN
            CALL BIMNXT ( LBCHCK )
            RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
            EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
         ELSE
            CALL NWINXT ( LBCHCK )
            IF ( BIMPCT .GT. RADTOT ) THEN
               NHOLE = NHOLE - 1
               ICYCL = ICYCL - 1
               CALL PHDSET ( IKPMX )
               IBRES = IBRES - IBAR (KPRIN)
               ICRES = ICRES - ICH  (KPRIN)
               BBRES = IBRES
               ZZRES = ICRES
               AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                ( BBRES, ZZRES)
               AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
               LTROUB = .FALSE.
               CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
               IF ( LTROUB ) THEN
                  KPNUCL (IKPMX) = 0
                  UMO2  = ERES**2 - PTRES2
                  UMO = SQRT (UMO2)
                  WRITE ( LUNOUT,* )' 0_P:UMO,AMNRES',UMO,AMNRES
                  IF ( KPRIN .EQ. 1 ) THEN
                     NPROT = NPROT + 1
                  ELSE
                     NNEUT = NNEUT + 1
                  END IF
                  LEXIT = .TRUE.
                  RETURN
               END IF
               EKNNUC = ENNUC (IKPMX) - AM (KPRIN)
               NP = NP + 1
               TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
               KPART (NP) = KPRIN
               PLR   (NP) = PNUCL  (IKPMX)
               CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
               CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
               CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
               WEI   (NP) = WEE
               KPNUCL (IKPMX) = 0
               IGREYP = IGREYP + ICH (KPRIN)
               IGREYN = IGREYN + 1 - ICH (KPRIN)
               PXINTR = PXINTR + PXNUCL (IKPMX)
               PYINTR = PYINTR + PYNUCL (IKPMX)
               PZINTR = PZINTR + PZNUCL (IKPMX)
               EINTR  = EINTR  + ENNUC  (IKPMX)
               IBINTR = IBINTR + IBAR   (KPART(NP))
               ICINTR = ICINTR + ICH    (KPART(NP))
               LEXIT  = .TRUE.
               RETURN
            END IF
            XSTNUC (IKPMX) = XIMPTR
            YSTNUC (IKPMX) = YIMPTR
            ZSTNUC (IKPMX) = ZIMPTR
            RSTNUC (IKPMX) = ABS (RIMPTR)
         END IF
         NHOLE = NHOLE - 1
         ICYCL = ICYCL - 1
         GO TO 1000
      END IF
      EKFNUC (NPNUC) = EKFPRO
      PXNUCL (NPNUC) = PCMSX + ETAX * PHELP
      PYNUCL (NPNUC) = PCMSY + ETAY * PHELP
      PZNUCL (NPNUC) = PCMSZ + ETAZ * PHELP
      PNUCL  (NPNUC) = SQRT ( PXNUCL (NPNUC)**2 + PYNUCL (NPNUC)**2
     &                      + PZNUCL (NPNUC)**2 )
      XSTNUC (NPNUC) = XIMPTR
      YSTNUC (NPNUC) = YIMPTR
      ZSTNUC (NPNUC) = ZIMPTR
      RSTNUC (NPNUC) = ABS (RIMPTR)
      RHNUCL (NPNUC) = RHOIMT
      NPNUC = NPNUC + 1
      KPNUCL (NPNUC) = KNUCIM
      KRFNUC (NPNUC) = KRFLIN + 1
      ETAPCM = - ETAPCM
      PHELP  = ECMSNU + ETAPCM / ( GAMCM + 1.D+00 )
      ENNUC  (NPNUC) = GAMCM * ECMSNU + ETAPCM
      IF ( ENNUC (NPNUC) - AM (KNUCIM) .LE. EKFIMP + DEFRNU ) THEN
         NPNUC  = NPNUC - NPLSIN
         LBCHCK = .FALSE.
         IF ( LBIMPC ) THEN
            CALL BIMNXT ( LBCHCK )
            RHOMEM = 0.5D+00 * ( RHOIMP + RHOIMT )
            EKFMEM = 0.5D+00 * ( EKFIMP + EKFPRO )
         ELSE
            CALL NWINXT ( LBCHCK )
            IF ( BIMPCT .GT. RADTOT ) THEN
               NHOLE = NHOLE - 1
               ICYCL = ICYCL - 1
               CALL PHDSET ( IKPMX )
               IBRES = IBRES - IBAR (KPRIN)
               ICRES = ICRES - ICH  (KPRIN)
               BBRES = IBRES
               ZZRES = ICRES
               AMMRES = BBRES * AMUAMU + 0.001D+00 * FKENER
     &                ( BBRES, ZZRES)
               AMNRES = AMMRES - ZZRES * AMELEC + ELBNDE ( ICRES )
               LTROUB = .FALSE.
               CALL UMOFIN ( IKPMX, BBRES, ZZRES, LTROUB )
               IF ( LTROUB ) THEN
                  KPNUCL (IKPMX) = 0
                  UMO2  = ERES**2 - PTRES2
                  UMO = SQRT (UMO2)
                  WRITE ( LUNOUT,* )' 0_T:UMO,AMNRES',UMO,AMNRES
                  IF ( KPRIN .EQ. 1 ) THEN
                     NPROT = NPROT + 1
                  ELSE
                     NNEUT = NNEUT + 1
                  END IF
                  LEXIT = .TRUE.
                  RETURN
               END IF
               EKNNUC = ENNUC (IKPMX) - AM (KPRIN)
               NP = NP + 1
               TKI   (NP) = ENNUC  (IKPMX) - AM (KPRIN)
               KPART (NP) = KPRIN
               PLR   (NP) = PNUCL  (IKPMX)
               CXR   (NP) = PXNUCL (IKPMX) / PLR (NP)
               CYR   (NP) = PYNUCL (IKPMX) / PLR (NP)
               CZR   (NP) = PZNUCL (IKPMX) / PLR (NP)
               WEI   (NP) = WEE
               KPNUCL (IKPMX) = 0
               IGREYP = IGREYP + ICH (KPRIN)
               IGREYN = IGREYN + 1 - ICH (KPRIN)
               PXINTR = PXINTR + PXNUCL (IKPMX)
               PYINTR = PYINTR + PYNUCL (IKPMX)
               PZINTR = PZINTR + PZNUCL (IKPMX)
               EINTR  = EINTR  + ENNUC  (IKPMX)
               IBINTR = IBINTR + IBAR   (KPART(NP))
               ICINTR = ICINTR + ICH    (KPART(NP))
               LEXIT = .TRUE.
               RETURN
            END IF
            XSTNUC (IKPMX) = XIMPTR
            YSTNUC (IKPMX) = YIMPTR
            ZSTNUC (IKPMX) = ZIMPTR
            RSTNUC (IKPMX) = ABS (RIMPTR)
         END IF
         NHOLE = NHOLE - 1
         ICYCL = ICYCL - 1
         GO TO 1000
      END IF
      EKFNUC (NPNUC) = EKFIMP
      PXNUCL (NPNUC) = -PCMSX + ETAX * PHELP
      PYNUCL (NPNUC) = -PCMSY + ETAY * PHELP
      PZNUCL (NPNUC) = -PCMSZ + ETAZ * PHELP
      PNUCL  (NPNUC) = SQRT ( PXNUCL (NPNUC)**2 + PYNUCL (NPNUC)**2
     &                      + PZNUCL (NPNUC)**2 )
      XSTNUC (NPNUC) = XIMPCT
      YSTNUC (NPNUC) = YIMPCT
      ZSTNUC (NPNUC) = ZIMPCT
      RSTNUC (NPNUC) = ABS (RIMPCT)
      RHNUCL (NPNUC) = RHOIMP
      POTOUT = ENNUC (NPNUC) - AM (KPRIN) + ENNUC (NPNUC-1) - AM(KNUCIM)
     &       - EKECON
      LBIMPC = .FALSE.
      LEXIT  = .FALSE.
      NUSCIN = NUSCIN + 1
      ISCTYP (NUSCIN) = KPRIN * 10 + KNUCIM
      NHLEXP = NHLEXP + 1
      HOLEXP (NHLEXP) = EKFIMP - EKFERM
      RHOEXP = RHOEXP + 0.5D+00 * ( RHOIMP + RHOIMT )
      EKFEXP = EKFEXP + 0.5D+00 * ( EKFIMP + EKFPRO )
      CALL NCLVFX
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
 3000 CONTINUE
      RETURN
*=== End of subroutine Nucnuc =========================================*
      END
