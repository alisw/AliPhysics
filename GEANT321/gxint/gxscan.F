*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXSCAN
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Scan Geometry control commands                            *
C.    *                                                                *
C.    *       Author:    R.Brun      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcscan.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcparm.inc"
      COMMON/QUEST/IQUEST(100)
      CHARACTER*20 CHOPT
      CHARACTER*32 CHPATL
      CHARACTER*4 NAME
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
      IF(CHPATL.EQ.'STURN')THEN
*
         IF(NPAR.LE.0) THEN
            IF(SCANFL) THEN
               CHOPT='ON'
            ELSE
               CHOPT='OFF'
            END IF
            WRITE(CHMAIL,10000) CHOPT
10000       FORMAT(' SCAN Parameters: SCAN mode is: ',A3)
            CALL GMAIL(0,0)
         ELSE
            CALL KUGETC(CHOPT,NCH)
            IF(CHOPT.EQ.'ON') THEN
               SCANFL=.TRUE.
            ELSE IF (CHOPT.EQ.'OFF') THEN
               SCANFL=.FALSE.
            ELSE IF (CHOPT.EQ.'INIT') THEN
               SCANFL=.TRUE.
               CALL GSCANI
            END IF
         END IF
      ELSEIF(CHPATL.EQ.'TETA')THEN
*
         IF(NPAR.LE.0) THEN
            WRITE(CHMAIL, 10100) NTETA, TETMIN, TETMAX, MODTET
10100       FORMAT(' SCAN Parameters: NTETA = ',I5,' MIN = ',G10.3,
     +             ' MAX = ',G10.3, ' MODE = ',I2)
            CALL GMAIL(0,0)
         ELSE
            CALL KUGETI(NTETA)
            CALL KUGETR(TETMIN)
            CALL KUGETR(TETMAX)
            IF(NSLMAX.EQ.0)NSLMAX=5
            CALL KUGETI(MODTET)
         END IF
      ELSEIF(CHPATL.EQ.'PHI')THEN
*
         IF(NPAR.LE.0) THEN
            WRITE(CHMAIL, 10200) NPHI, PHIMIN, PHIMAX
10200       FORMAT(' SCAN Parameters: NPHI  = ',I5,' MIN = ',G10.3,
     +             ' MAX = ',G10.3)
            CALL GMAIL(0,0)
         ELSE
            CALL KUGETI(NPHI)
            CALL KUGETR(PHIMIN)
            CALL KUGETR(PHIMAX)
            IPHI1=1
            IPHIL=NPHI
            IF(NSLMAX.EQ.0)NSLMAX=5
         END IF
      ELSEIF(CHPATL.EQ.'SLIST')THEN
*
         IF(NPAR.LE.0) THEN
            NPADON=MIN(NSLIST,8)
            NPAWRI=NPADON
            WRITE(CHMAIL,10300) NSLIST, (ISLIST(J), J=1, NPAWRI)
10300       FORMAT(' SCAN Parameters: ',I3,' Scan volumes :',8(1X,A4))
            CALL GMAIL(0,0)
   10       IF(NSLIST.GT.NPADON) THEN
               NPAWRI=MIN(NSLIST-NPADON,15)
               WRITE(CHMAIL,10400)(ISLIST(J),J=NPADON+1,NPADON+NPAWRI)
10400          FORMAT((1X,15(1X,A4)))
               CALL GMAIL(0,0)
               NPADON=NPADON+NPAWRI
               GO TO 10
            ENDIF
         ELSE
            CALL VBLANK(ISLIST, MSLIST)
            IF(NPAR.GT.MSLIST) THEN
               WRITE(CHMAIL,10500) MSLIST
10500          FORMAT(' Warning! Only first ',I3,' scan volumes ',
     +                'accepted')
               CALL GMAIL(0,0)
            ENDIF
            NSLIST=0
            NVOL=IQ(JVOLUM-1)
            DO 20 I=1,MIN(NPAR,MSLIST)
               CALL KUGETC(CHOPT,NCH)
               IF(I.EQ.1) THEN
                  IF(NPAR.EQ.1.AND.CHOPT.EQ.'.') THEN
                     CALL MZDROP(IXCONS,LSCAN,' ')
                     GO TO 999
                  END IF
               END IF
               NCH=MIN(4,NCH)
               NSLIST=NSLIST+1
               CALL UCTOH(CHOPT,ISLIST(NSLIST),4,NCH)
               JVOL=IUCOMP(ISLIST(NSLIST),IQ(JVOLUM+1),NVOL)
               IF(JVOL.LE.0) THEN
                  WRITE(CHMAIL,10600) ISLIST(NSLIST)
10600             FORMAT(' Warning: volume ',A4,' does not exist;',
     +                   ' skipped')
                  CALL GMAIL(0,0)
                  NSLIST=NSLIST-1
               END IF
   20       CONTINUE
            IF(NSLIST.LE.0) THEN
               WRITE(CHMAIL,10700)
10700          FORMAT(' Warning! No valid volume defined')
               CALL GMAIL(0,0)
            ENDIF
         ENDIF
      ELSEIF(CHPATL.EQ.'VERTEX')THEN
*
         IF(NSLMAX.EQ.0)NSLMAX=5
         CALL KUGETR(VSCAN(1))
         CALL KUGETR(VSCAN(2))
         CALL KUGETR(VSCAN(3))
      ELSEIF(CHPATL.EQ.'PCUTS')THEN
*
         IF(NPAR.LE.0) THEN
            WRITE(CHMAIL,10800) IPARAM
10800       FORMAT(' Parametrization flag = ',I2,
     +             ' Parametrization cuts:')
            CALL GMAIL(0,0)
            WRITE(CHMAIL,10900)
10900       FORMAT('         Gamma     Electrons   Ch. Hadrons',
     +             '  Neu. Hadrons         Muons')
            CALL GMAIL(0,0)
            WRITE(CHMAIL,11000) PACUTS
11000       FORMAT(5(1X,G13.4))
            CALL GMAIL(0,0)
         ELSE
            CALL KUGETI(IPARAM)
            DO 30 JPACUT=1, 5
               CALL KUGETR(PACUTS(JPACUT))
   30       CONTINUE
         ENDIF
      ELSEIF(CHPATL.EQ.'SFACTORS')THEN
*
         CALL KUGETR(FACTX0)
         CALL KUGETR(FACTL)
         CALL KUGETR(FACTR)
      ELSEIF(CHPATL.EQ.'LSCAN')THEN
*
         CALL KUGETI(IDPHI)
         CALL KUGETC(NAME,NCH)
         CALL KUGETC(CHOPT,NCH)
         CALL GXSCAL(IDPHI,NAME,CHOPT)
         IF(INDEX(CHOPT,'P').NE.0)THEN
            IF(IDPHI.NE.0)CALL HPLEGO(IDPHI,30.,30.)
         ENDIF
      ELSEIF(CHPATL.EQ.'HSCAN')THEN
*
         CALL KUGETI(IDPHI)
         CALL KUGETC(NAME,NCH)
         CALL KUGETC(CHOPT,NCH)
         CALL GXSCAH(IDPHI,NAME,CHOPT)
         IF(INDEX(CHOPT,'P').NE.0)THEN
            IF(IDPHI.NE.0)CALL HPLOT(IDPHI,' ',' ',0)
         ENDIF
      ENDIF
*
  999 END
