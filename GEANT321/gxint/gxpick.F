*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXPICK
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Geometry commands                                         *
C.    *                                                                *
C.    * Point to volume just drawn to pick up medium name, volume name *
C.    * etc. The first point points to the volume, the second point    *
C.    * gives the position of the character string which contains this *
C.    * information.                                                   *
C.    *                                                                *
C.    *       Authors:   S.Egli      **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gconsp.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcbank.inc"
#include "geant321/gcsets.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcmate.inc"
#include "geant321/gcunit.inc"
*
* this COMMON filled in routine GDRAW !
*
      COMMON/GCVHLP/NVLAST
*
      DIMENSION XC(3),XYZ(3),XINVMA(3,3),VL(3),VM(3),VN(3)
      DIMENSION XX(10),YY(10),RHELP(3)
*
      CHARACTER*4 NAMV,CHIDTY
      CHARACTER*20 NAMM
      LOGICAL BTEST
*
* determine inverse matrix xinvma for current view parameters
*
      PH = ABS(MOD(GPHI,360.))
      THET = ABS(MOD(GTHETA,360.))
      IF(THET.LE.180.)GO TO 10
      PH = PH + 180.
      THET = 360. - THET
*
   10 ST = SIN(THET * DEGRAD)
      CT = COS(THET * DEGRAD)
      SP = SIN(PH * DEGRAD)
      CP = COS(PH * DEGRAD)
*
*             VN is new nu axis
*
      VN(1) = ST * CP
      VN(2) = ST * SP
      VN(3) = CT
*
      IF(ABS(VN(2)).GT.0.99999) THEN
*
*             Special case when observer line of sight is along mu:
*             in this case one chooses arbitrarily the vertical axis of
*             plane of projection as the lambda axis and the horizontal
*             as the nu axis
*
         VL(1) = 0.
         VL(2) = 0.
         VL(3) = 1.
         VM(1) = 1.
         VM(2) = 0.
         VM(3) = 0.
      ELSE
*
         VM(1) = 0.
         VM(2) = 1.
         VM(3) = 0.
*
*             Define new lambda axis
*
         CALL CROSS(VM,VN,VL)
         CALL VUNIT(VL,VL,3)
*
*             Define new mu axis
*
         CALL CROSS(VN,VL,VM)
      ENDIF
*
*   now invert matrix defined by VL,VM,VN -> XINVMA
*
      DO 20 I=1,3
         XINVMA(1,I)=VL(I)
         XINVMA(2,I)=VM(I)
         XINVMA(3,I)=VN(I)
   20 CONTINUE
*
      CALL RINV(3,XINVMA,3,RHELP,IFAIL)
      IF(IFAIL.NE.0)THEN
         WRITE(CHMAIL,10100) IFAIL
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
*
*   perspective projection ?
*
      CALL UCTOH('PERS',IPERS,4,4)
      IF(IPRJ.EQ.IPERS)THEN
         WRITE(CHMAIL,10200)
         CALL GMAIL(0,0)
         GOTO 999
      ENDIF
*
* pick up two points in user coordinates:
*
   30 CALL IRQLC(1,1,ISTAT,NT,U0,V0)
      IF(ISTAT.EQ.0)GOTO 999
      CALL IRQLC(1,1,ISTAT,NT,U1,V1)
      IF(ISTAT.EQ.0)GOTO 999
*
* transform (u0,v0) to coordinates in MARS system:
*  (inverse operation of what is done in routine GDFR3D)
*
*   take zoom parameters into account:
*
      U01=(U0-GZUB-GZUC)/GZUA
      V01=(V0-GZVB-GZVC)/GZVA
*
*   rotate and shift back
*
      UU=+COSPSI*(U01-GU0)+SINPSI*(V01-GV0)
      VV=-SINPSI*(U01-GU0)+COSPSI*(V01-GV0)
      XYZ(1)=UU/GSCU
      XYZ(2)=VV/GSCV
      XYZ(3)=DCUT
*
*   apply xinvma
*
      XC(1)=XINVMA(1,1)*XYZ(1)+XINVMA(1,2)*XYZ(2)+XINVMA(1,3)*XYZ(3)
      XC(2)=XINVMA(2,1)*XYZ(1)+XINVMA(2,2)*XYZ(2)+XINVMA(2,3)*XYZ(3)
      XC(3)=XINVMA(3,1)*XYZ(1)+XINVMA(3,2)*XYZ(2)+XINVMA(3,3)*XYZ(3)
*
*   build up GCVOLU structure with last drawn volume as
*   top of tree
*
      NLEV=1
      LNUM=0
      CALL GLVOLU(NLEV,NVLAST,LNUM,IER)
*
*  determine medium
*
      NUMED=0
      CALL GMEDIA(XC,NUMED)
*
      IF(NUMED.EQ.0)THEN
         WRITE(CHMAIL,10300)
         CALL GMAIL(0,0)
         GOTO 30
      ENDIF
      JTM = LQ(JTMED- NUMED)
      DO 40 I=1,5
   40 NATMED(I)=IQ(JTM+I)
      NMAT   = Q(JTM + 6)
      ISVOL  = Q(JTM + 7)
      IFIELD = Q(JTM + 8)
      FIELDM = Q(JTM + 9)
      TMAXFD = Q(JTM + 10)
      STEMAX = Q(JTM + 11)
      DEEMAX = Q(JTM + 12)
      EPSIL  = Q(JTM + 13)
      STMIN  = Q(JTM + 14)
      CALL UHTOC(NAMES(NLEVEL),4,NAMV,4)
      CALL UHTOC(NATMED,4,NAMM,20)
      DO 50 I=1,20
         IF(NAMM(I:I).EQ.'$')NAMM(I:I)=' '
   50 CONTINUE
      DO 60 I=20,1,-1
         IF(NAMM(I:I).NE.' ')GOTO 70
   60 CONTINUE
   70 NJLAST=I
*
*  determine detector idtype
*
      IF(JSET.GT.0)CALL GFINDS
*
* draw pointer and write volume name,medium
*
      XX(1)=U0
      YY(1)=V0
      XX(2)=U1
      YY(2)=V1
      CALL IPL(2,XX,YY)
*
*  determine text alignment
*
      PHI=ATAN2(V1-V0,U1-U0)*RADDEG
      IF(ABS(PHI).LT.90.)THEN
         IHOR=1
      ELSE
         IHOR=3
      ENDIF
      IF(PHI.GT.0.)THEN
         IVER=5
      ELSE
         IVER=1
      ENDIF
      CALL ISTXAL(IHOR,IVER)
*
      IF(.NOT.BTEST(IQ(LQ(JVOLUM-LVOLUM(1))),4))THEN
         CALL ITX(U1,V1,NAMV//'"j# '//NAMM(1:NJLAST))
      ELSE
         WRITE(CHIDTY,10000)IDTYPE
         CALL ITX(U1,V1,NAMV//'"j# '//NAMM(1:NJLAST)//','//CHIDTY)
      ENDIF
      GOTO 30
10000 FORMAT(I4)
10100 FORMAT(' GXPICK: Matrix inversion failed with ',I3,
     +       '; abandoning')
10200 FORMAT(' GXPICK: perspective projection can not be handled')
10300 FORMAT(' GXPICK: point is outside volume')
  999 END
