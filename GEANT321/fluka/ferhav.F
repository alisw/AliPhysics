*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 02/07/94  17.46.10  by  S.Giani
*-- Author :
      SUBROUTINE FERHAV ( KP, EPROJ, PPROJ, TXX, TYY, TZZ )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
#include "geant321/balanc.inc"
#include "geant321/finlsp.inc"
#include "geant321/hadflg.inc"
#include "geant321/nucdat.inc"
#include "geant321/qquark.inc"
#include "geant321/part3.inc"
#include "geant321/resnuc.inc"
      COMMON /FKABLT/ AM(110), GAA(110), TAU(110), ICH(110), IBAR(110),
     &                K1(110), K2(110)
      COMMON / FKNUCF / DELEFT, EKRECL, V0EXTR, ITTA, ITJ, LVMASS
      LOGICAL LVMASS, LSMPAN
      COMMON / FKEVNT / LNUCRI, LHADRI
      LOGICAL LNUCRI, LHADRI
      REAL FRNDM(3)
      SAVE ONEDUM, ZERDUM
      DATA ONEDUM / 1.D+00 /
      DATA ZERDUM / 0.D+00 /
      AMPROJ = AM (KP)
      AMTAR  = AM (ITTA)
      ECHCK  = EPROJ + EFRM - V0WELL (ITJ) - EKRECL - EBNDNG (ITJ)
      PXCHCK = PXFRM + PPROJ * TXX
      PYCHCK = PYFRM + PPROJ * TYY
      PZCHCK = PZFRM + PPROJ * TZZ
      UMIN2  = ( AMPROJ + AMTAR )**2
      P2CHCK = PXCHCK**2 + PYCHCK**2 + PZCHCK**2
      UMO2   = ECHCK**2  - P2CHCK
      IF ( ABS ( UMO2 - UMIN2 ) .GT. TWOTWO * ANGLGB * UMO2 ) THEN
         EPROJX = HLFHLF * ( UMO2 - AMPROJ**2 - AMTAR**2 ) / AMTAR
         IF ( EPROJX .LT. AMPROJ ) THEN
            WRITE (LUNERR,*)' Ferhav: trouble with pseudo-masses!!',
     &      EPROJX,AMPROJ,LVMASS
            EPROJX = AMPROJ
            PPROJX = ZERZER
            LRESMP = .TRUE.
            RETURN
         ELSE
            PPROJX = SQRT ( ( EPROJX - AMPROJ ) * ( EPROJX + AMPROJ ) )
         END IF
         ETOTX  = EPROJX + AMTAR
         AMTRMX = HLFHLF * ( UMO2 - AMPROJ**2 - AMTAR**2 ) / AMPROJ
         AMSQMX = AMTRMX**2
         PTOSCA = PXCHCK * TXX + PYCHCK * TYY + PZCHCK * TZZ
         PXTART = PXCHCK - PTOSCA * TXX
         PYTART = PYCHCK - PTOSCA * TYY
         PZTART = PZCHCK - PTOSCA * TZZ
         PTRASQ = PXTART**2 + PYTART**2 + PZTART**2
         AMTRSQ = AMTAR**2  + PTRASQ
         IF ( AMTRSQ .GT. AMSQMX ) THEN
            PPCHCK = SQRT (P2CHCK)
            PTOOLD = PTOSCA
            PTRASQ = ( AMTRMX - AMTAR ) * ( AMTRMX + AMTAR )
******            AMTRSQ = AMTRMX + PTRASQ
******            PTOSCA = SIGN (ONEONE,PTOOLD) * SQRT ( P2CHCK - PTRASQ )
            AMTRSQ = AMSQMX
            PTOSCA = SQRT ( P2CHCK - PTRASQ )
            ALPTUU = ( PTOSCA * SQRT ( ( PPCHCK - PTOOLD ) * ( PPCHCK
     &             + PTOOLD ) / ( PPCHCK - PTOSCA ) / ( PPCHCK + PTOSCA
     &             ) ) - PTOOLD ) / PPCHCK
            FNORM  = SQRT ( ONEONE + ALPTUU**2 + TWOTWO * ALPTUU
     &             * PTOOLD / PPCHCK )
            ALPTUU = ALPTUU / PPCHCK
            TXXX   = ( TXX + ALPTUU * PXCHCK ) / FNORM
            TYYY   = ( TYY + ALPTUU * PYCHCK ) / FNORM
            TZZZ   = ( TZZ + ALPTUU * PZCHCK ) / FNORM
            UMOTR2 = UMO2 + PTRASQ
            UMOTR  = SQRT (UMOTR2)
            AMTRAN = AMTRMX
            PPARSQ = PTOSCA**2
            PPARTT = PTOSCA
            GAMCMS = ECHCK  / UMOTR
            ETACMS = PPARTT / UMOTR
            EPRCMS = AMPROJ
            PPRCMS = ZERZER
         ELSE
            TXXX = TXX
            TYYY = TYY
            TZZZ = TZZ
            PTOOLD = PTOSCA
            UMOTR2 = UMO2 + PTRASQ
            UMOTR  = SQRT (UMOTR2)
            PPARTT = PTOSCA
            GAMCMS = ECHCK  / UMOTR
            ETACMS = PPARTT / UMOTR
            EPRCMS = HLFHLF * ( UMOTR2 + AMPROJ**2 - AMTRSQ ) / UMOTR
            PPRCMS = SQRT ( ( EPRCMS - AMPROJ ) * ( EPRCMS + AMPROJ ) )
         END IF
         EPRLAB = GAMCMS * EPRCMS + ETACMS * PPRCMS
         ETRLAB = ECHCK  - EPRLAB
         PPRLAB = SQRT ( ( EPRLAB - AMPROJ ) * ( EPRLAB + AMPROJ ) )
         PXTARG = PXCHCK - PPRLAB * TXXX
         PYTARG = PYCHCK - PPRLAB * TYYY
         PZTARG = PZCHCK - PPRLAB * TZZZ
         GAM    = ETRLAB / AMTAR
         BGX    = PXTARG / AMTAR
         BGY    = PYTARG / AMTAR
         BGZ    = PZTARG / AMTAR
         PPHELP = ( BGX * TXXX + BGY * TYYY + BGZ * TZZZ ) * PPRLAB
         ETAPCM = EPRLAB - PPHELP / ( GAM + ONEONE )
         PXPROJ = PPRLAB * TXXX - BGX * ETAPCM
         PYPROJ = PPRLAB * TYYY - BGY * ETAPCM
         PZPROJ = PPRLAB * TZZZ - BGZ * ETAPCM
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
      ELSE
         UMO2   = UMIN2
         EPROJX = AMPROJ
         PPROJX = ZERZER
         ETOTX  = AMPROJ + AMTAR
         LSMPAN = .FALSE.
         CALL POLI   ( COSTH0, SINTH0 )
         CALL SFECFE ( SINPH0, COSPH0 )
         UUOLD  = SINTH0 * COSPH0
         VVOLD  = SINTH0 * SINPH0
         WWOLD  = COSTH0
         AAFACT = ECHCK  + ETOTX
         BBFACT = PPROJX - PZCHCK
         DDENOM = ETOTX * AAFACT - PPROJX * BBFACT
         GAM = ( ECHCK * AAFACT + PPROJX * BBFACT ) / DDENOM
         BGZ = - BBFACT * AAFACT / DDENOM
         BGX = PXCHCK * ( GAM + ONEONE ) / AAFACT
         BGY = PYCHCK * ( GAM + ONEONE ) / AAFACT
      END IF
      PLABS  = PPROJX
      ELABS  = EPROJX
      IF ( PLABS .LT. 1.D-04 ) THEN
         WRITE (LUNERR,*)' Ferhav: kp,plabs,elabs,pprox,y,z,pfrmix,y,z'
     &   ,KP,PLABS,ELABS,PPROJ*TXX,PPROJ*TYY,PPROJ*TZZ,PXFRM,PYFRM,
     &    PZFRM
         WRITE (LUNERR,*)'   Lvmass,Am(kp),Eproj:',LVMASS,AM(KP),EPROJ
      ELSE IF ( PLABS .GT. 1.D+01 ) THEN
         WRITE (LUNERR,*)' Ferhav: kp,plabs,elabs,pprox,y,z,pfrmix,y,z'
     &   ,KP,PLABS,ELABS,PPROJ*TXX,PPROJ*TYY,PPROJ*TZZ,PXFRM,PYFRM,
     &    PZFRM
         WRITE (LUNERR,*)'   Lvmass,Am(kp),Eproj:',LVMASS,AM(KP),EPROJ
      END IF
      ISSU = 0
      DO 100 IQ = 1,3
         ISSU = ISSU + MQUARK (IQ,KP) / 3
  100 CONTINUE
      IF ( LVMASS ) THEN
         LHADRI = .TRUE.
         CALL HADRIN ( KP, PLABS, ELABS, ZERDUM, ZERDUM, ONEDUM, ITTA )
         IOLDHD = 0
      ELSE IF ( PLABS .GT. 7.D+00 ) THEN
         LHADRI = .FALSE.
         CALL HINHEV ( KP, PLABS, ELABS, ITTA )
      ELSE
         LHADRI = .TRUE.
         CALL HADRIV ( KP, PLABS, ELABS, ZERDUM, ZERDUM, ONEDUM, ITTA )
      END IF
      DO 2000 I=1,IR
         ECMS  = ELR (I)
         PCMSX = PLR (I) * CXR (I)
         PCMSY = PLR (I) * CYR (I)
         PCMSZ = PLR (I) * CZR (I)
         IF ( LSMPAN ) THEN
            PCMSX = PLR (I) * CXR (I)
            PCMSY = PLR (I) * CYR (I)
            PCMSZ = WWOLD * PLR (I) * CZR (I)
         ELSE
            PLRX = CXR (I) * COSPH0 * WWOLD - CYR (I) * SINPH0
     &           + CZR (I) * UUOLD
            PLRY = CXR (I) * SINPH0 * WWOLD + CYR (I) * COSPH0
     &           + CZR (I) * VVOLD
            PLRZ = - CXR (I) * SINTH0 + CZR (I) * WWOLD
            PCMSX = PLRX * PLR (I)
            PCMSY = PLRY * PLR (I)
            PCMSZ = PLRZ * PLR (I)
         END IF
         CALL ALTRA ( GAM, BGX, BGY, BGZ, PCMSX, PCMSY, PCMSZ,
     &                ECMS, PLR (I), PLRX, PLRY, PLRZ, ELR (I) )
         CXR (I) = PLRX / PLR (I)
         CYR (I) = PLRY / PLR (I)
         CZR (I) = PLRZ / PLR (I)
         DO 200 IQ = 1,3
            ISSU = ISSU - MQUARK (IQ,KPTOIP(ITR(I))) / 3
  200    CONTINUE
2000  CONTINUE
      IF ( ISSU .NE. 0 ) THEN
         WRITE (LUNOUT,*)' *** Strangeness non conservation in Hadriv',
     &                     ISSU,KP,ITTA,' ***'
         WRITE (LUNERR,*)' *** Strangeness non conservation in Hadriv',
     &                     ISSU,KP,ITTA,' ***'
         LRESMP = .TRUE.
      END IF
      V0WELL (ITJ) = V0WELL (ITJ) - V0EXTR
      RETURN
      ENTRY FERSET
      FERM = PFRMMX (ITJ)
      CALL GRNDM(FRNDM,3)
      P2 = MAX ( FRNDM (1), FRNDM (2), FRNDM (3) )
      IF ( IBTAR .LE. 1 ) THEN
         FERM = ZERZER
      END IF
      P2=FERM*P2
      P2SQ   = P2 * P2
      IATEMP = IBTAR - 1
      ATEMP  = DBLE ( IBTAR ) - ONEONE
      IF ( ITJ .EQ. 1 ) THEN
         IZTEMP = ICHTAR - 1
      ELSE
         IZTEMP = ICHTAR
      END IF
      ZTEMP =  DBLE ( IZTEMP )
      DELCTR = ( DBLE (ICHTAR) - ZTEMP ) * AMELEC
      DELEFT = AMMTAR - AMNTAR - DELCTR
      AMMRES = AMUAMU * ATEMP + 1.D-03 * FKENER ( ATEMP, ZTEMP )
      AMNRES = AMMRES - ZTEMP * AMELEC + ELBNDE ( IZTEMP )
      AMTMSQ = AMMRES * AMMRES
      EKRECL = SQRT ( AMTMSQ + P2SQ ) - AMMRES
      CALL POLI ( POLC, POLS )
      CALL COSI ( SFE,  CFE )
      PXFRM = CFE * POLS * P2
      PYFRM = SFE * POLS * P2
      PZFRM = POLC * P2
      EFRM  = SQRT ( AMNUSQ (ITJ) + P2SQ )
      EKFER  = EFRM  - AMNUCL (ITJ)
      TVEUZ  = V0WELL (ITJ) - EFRM + EBNDNG (ITJ) + AMMTAR - AMMRES
     &       - DELCTR
      IF ( TVEUZ .LT. ZERZER ) THEN
         V0EXTR = - TVEUZ + TENTEN * TVEPSI
         TVEUZ  = TVEUZ  + V0EXTR
         V0WELL (ITJ) = V0WELL (ITJ) + V0EXTR
      ELSE
         V0EXTR = ZERZER
      END IF
      RETURN
      END
