*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:25  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDRAWT(X,Y,CHTEXT,SIZE,ANGLE,LWIDTH,IOPT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Draw text with software characters                       *
C.    *                                                                *
C.    *       In ITXATT (common GCDRAW) there is                       *
C.    *       a bit mask for the text attributes :                     *
C.    *                                                                *
C.    *       Bit  1- 3 = Text line width (LWIDTH)                     *
C.    *       Bit  4- 7 = Text color                                   *
C.    *       Bit  8-11 = Text font                                    *
C.    *       Bit 12-24 = Reserved for future use                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDAHIT, GDARR, GDHEAD,      *
C.    *                   GDPART, GDPRTR, GDSCAL, GDSHOW, GDSPEC       *
C.    *       Authors : P.Zanarini ; S.Giani  *********                *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gcflag.inc"
      CHARACTER*(*) CHTEXT
      CHARACTER*1  CHOPT
      DIMENSION AITXFO(2)
C.
C.    ------------------------------------------------------------------
C.
      CALL IGQ('TXCI',AITXCO)
      CALL IGQ('TXFP',AITXFO)
      ALWIDT=LWIDTH
      NCHAR=LNBLNK(CHTEXT)
      IF (IDVIEW.EQ.0.OR.IDVIEW.EQ.-175) GO TO 40
C
C             Store text on view bank IDVIEW
C
      JV=LQ(JDRAW-IDVIEW)
      ITU=ITU+1
C
   10 IF (ITU.LE.MAXTU) GO TO 20
C
C             Push text unit bank
C
      JV = LQ(JV-3)
      CALL MZPUSH(IXCONS,JV,0,MORTU,'I')
      IF(IEOTRI.NE.0)GO TO 50
      JV=LQ(JDRAW-IDVIEW)
      MAXTU=MAXTU+MORTU
      GO TO 10
C
   20 IF ((ITS+6+NCHAR).LE.MAXTS) GO TO 30
C
C             Push text segment bank
C
      JV = LQ(JV-6)
      CALL MZPUSH(IXCONS,JV,0,MORTS,'I')
      IF(IEOTRI.NE.0)GO TO 50
      JV=LQ(JDRAW-IDVIEW)
      MAXTS=MAXTS+MORTS
      GO TO 20
C
   30 JVTS=LQ(JV-6)+ITS+1
      Q(JVTS)=X
      Q(JVTS+1)=Y
      Q(JVTS+2)=SIZE
      Q(JVTS+3)=ANGLE
C
C             Bit 1- 3 = Text line width (LWIDTH)
C             Bit 4-24 = Text attribute
C
      ISUM=0
      ITXCOL=AITXCO
      CALL MVBITS(ITXCOL,0,4,ITXATT,3)
      IF(AITXFO(1).EQ.-6..AND.AITXFO(2).EQ.0.)THEN
         ITXFON=6
      ELSEIF(AITXFO(1).EQ.-6..AND.AITXFO(2).EQ.1.)THEN
         ITXFON=7
      ENDIF
      CALL MVBITS(ITXFON,0,4,ITXATT,7)
      CALL MVBITS(ITXATT,0,24,ISUM,0)
      CALL MVBITS(LWIDTH,0,3,ISUM,0)
      Q(JVTS+4)=ISUM
C
      Q(JVTS+5)=IOPT
      Q(JVTS+6)=NCHAR
      CALL UCTOH(CHTEXT,IQ(JVTS+7),4,NCHAR)
C
      Q(LQ(JV-3)+ITU)=ITS+1
      ITS=ITS+7+NCHAR
      GO TO 999
C
C             Draw text on screen
C
   40 CONTINUE
      IF (IOPT.EQ.-1) THEN
         CHOPT='L'
         ALFLAG=10.
      ELSE IF (IOPT.EQ.0) THEN
         CHOPT='C'
         ALFLAG=20.
      ELSE IF (IOPT.EQ.1) THEN
         CHOPT='R'
         ALFLAG=30.
      ENDIF
      IF(IDVIEW.EQ.-175)THEN
         ITXCOL=IBITS(ITXATT,3,4)
         AITXCO=ITXCOL
         ITXFON=IBITS(ITXATT,7,4)
         IF(ITXFON.EQ.6)THEN
            ITXFON=-60
         ELSEIF(ITXFON.EQ.7)THEN
            ITXFON=-61
         ENDIF
         AITXFP=ITXFON
         CALL IGSET('TXFP',AITXFP)
      ENDIF
      IF(ALWIDT.GT.1.)ALWIDT=ALWIDT*3.
***      CALL IGTEXT(X,Y,CHTEXT,SIZE,ANGLE,CHOPT)
      CALL IGSET('TXAL',ALFLAG)
      CALL IGSET('TANG',ANGLE)
      CALL IGSET('TXCI',AITXCO)
      CALL IGSET('CHHE',SIZE)
      CALL IGSET('LWID',ALWIDT)
      CALL ITX(X,Y,CHTEXT)
      GO TO 999
C
   50 WRITE (CHMAIL,10000)
      CALL GMAIL(0,0)
C
10000 FORMAT (' ERROR IN GDRAWT: MEMORY OVERFLOW IN PUSHING A BANK')
  999 END
