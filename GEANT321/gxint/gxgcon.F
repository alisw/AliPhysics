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
      SUBROUTINE GXGCON
C.
C.    ******************************************************************
C.    *                                                                *
C.    *      Graphics control commands                                 *
C.    *                                                                *
C.    *       Authors:   R.Brun      **********                        *
C.    *                  P.Zanarini  **********                        *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcxlun.inc"
*KEEP,HIMETA
      COMMON /HIMETA/ IDMETA,XMFACT,YMFACT,TEKACT,METACT,FILOPN
      LOGICAL TEKACT,METACT,FILOPN
      COMMON /HIMETC/ CHMETA
      CHARACTER*24 CHMETA
 
      COMMON/CWK/IWK
      COMMON/DEFPAR/IATDEF,IVADEF
      CHARACTER*4 CHNAME,CHOPT,CHPAR,IATDEF
      CHARACTER*32 CHPATL
      LOGICAL IGM100
C.
C.    ------------------------------------------------------------------
C.
      CALL KUPATL(CHPATL,NPAR)
*
C             METAFILE
C
      IF(CHPATL.EQ.'METAFILE')THEN
         CALL KUGETI(I)
         LUN=IABS(I)
         CALL KUGETI(METAFL)
         CALL KUGETS(CHMETA,NCH)
         IF(IGM100(METAFL))THEN
            IF(LUN.LT.100)THEN
               IF(I.LT.0)THEN
                  I=I-100
               ELSEIF(I.GT.0)THEN
                  I=I+100
               ENDIF
            ELSEIF(LUN.GT.100)THEN
               LUN=LUN-100
            ENDIF
         ENDIF
         IF(LUN.NE.0)THEN
            CALL GXLUNF(LUN,2,IFREE)
            IF(IFREE.EQ.0)GO TO 99
         ENDIF
         LUNIT(LUN)=8
         IF(IWK.LE.0)THEN
            IF(I.GT.0)THEN
               IWK=-2
            ELSEIF(I.LT.0)THEN
               IWK=-1
            ELSE
               IWK=0
            ENDIF
         ENDIF
         CALL IGMETA(I,METAFL)
*
      ELSEIF (CHPATL.EQ.'DOPEN') THEN
         CALL KUGETI(IVIEW)
         CALL GDOPEN(IVIEW)
*
      ELSEIF (CHPATL.EQ.'DSHOW') THEN
         CALL KUGETI(IVIEW)
         IF(IVIEW.EQ.1000)THEN
           IVIEW=MYISEL
         ELSE
           MYISEL=IVIEW
         ENDIF
         CALL GDSHOW(IVIEW)
*
      ELSEIF (CHPATL.EQ.'DELETE') THEN
         CALL KUGETI(IVIEW)
         CALL GDELET(IVIEW)
*
      ELSEIF (CHPATL.EQ.'DCLOSE') THEN
         CALL GDCLOS
*
      ELSEIF (CHPATL.EQ.'CHANGEWK') THEN
         CALL GCHNWK
*
      ELSEIF (CHPATL.EQ.'RESETWK') THEN
         CALL GRESWK
*
      ELSEIF (CHPATL.EQ.'SSETAT') THEN
         CALL KUGETC(CHOPT,NCH)
         IATDEF=CHOPT
*
      ELSEIF (CHPATL.EQ.'SSETVA') THEN
         CALL KUGETR(RVAL)
         IVADEF=RVAL
*
      ELSEIF (CHPATL.EQ.'SATT') THEN
         CALL KUGETC(CHNAME,NCH)
         CALL KUGETC(CHOPT,NCH)
         CALL KUGETI(IVAL)
         CALL GSATT(CHNAME,CHOPT,IVAL)
*
      ELSEIF (CHPATL.EQ.'SCALE') THEN
         CALL KUGETR(GSCU)
         CALL KUGETR(GSCV)
*
      ELSEIF (CHPATL.EQ.'COLOR') THEN
         CALL KUGETI(ICOL)
         ICOL=-IABS(ICOL)
         CALL GDCOL(ICOL)
*
      ELSEIF (CHPATL.EQ.'LWID') THEN
         CALL KUGETI(LW)
         LW=-IABS(LW)
         CALL GDLW(LW)
*
      ELSEIF (CHPATL.EQ.'NEXT') THEN
         CALL HPLI
         CALL ISELNT(1)
         CALL IGRNG(20.,20.)
*
      ELSEIF (CHPATL.EQ.'DOPT') THEN
         CALL KUGETC(CHOPT,NCH)
         CALL KUGETC(CHPAR,NCH)
         CALL GDOPT(CHOPT,CHPAR)
*
      ELSEIF (CHPATL.EQ.'SIZE') THEN
         CALL KUGETR(PLT1)
         PLT2=PLT1
         CALL KUGETR(PLT2)
         PLTRNX=PLT1
         PLTRNY=PLT2
         CALL IGRNG(PLT1,PLT2)
         CALL HPLSET('XSIZ',PLT1)
         CALL HPLSET('YSIZ',PLT2)
*
      ELSEIF (CHPATL.EQ.'SPERS') THEN
         CALL KUGETR(DPERS)
         DPERS=ABS(DPERS)
*
      ELSEIF (CHPATL.EQ.'MAP_COLOR') THEN
         IADD=0
         IF (NPAR.GT.0) THEN
   10       CALL KUGETI(IADD)
            IF (IADD.LT.0.OR.IADD.GT.16) GO TO 10
         ENDIF
         IF (NPAR.GT.1) THEN
   20       CALL KUGETI(IVAL)
            IF (IVAL.LT.1.OR.IVAL.GT.16) GO TO 20
         ENDIF
         IF (IADD.EQ.0) THEN
            IAD1=1
            IAD2=16
         ELSE
            IAD1=IADD
            IAD2=IADD
         ENDIF
         DO 30 I=IAD1,IAD2
            IF (NPAR.GT.1) LOOKTB(I)=IVAL
            WRITE (CHMAIL,'('' COLOR('',I2,'') = '',I2)') I,LOOKTB(I)
            CALL GMAIL(0,0)
   30    CONTINUE
      ENDIF
*
  99  END
