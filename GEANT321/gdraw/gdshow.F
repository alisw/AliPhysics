*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:28  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDSHOW(IVIEW)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Show view bank IVIEW on the screen                       *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>                              *
C.    *       Authors : R.Brun, P.Zanarini, S.Giani   *********        *
C.    *                                                                *
C.    ******************************************************************
C.
C
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gcunit.inc"
#include "geant321/gconst.inc"
      CHARACTER*80 CHTEXT
C.
C.    ------------------------------------------------------------------
C.
      IF (IVIEW.LE.0) GO TO 70
      IF (IVIEW.GT.NKVIEW) GO TO 70
      IF (LQ(JDRAW-IVIEW).EQ.0) GO TO 70
C
      JV=LQ(JDRAW-IVIEW)
      GTHETA=Q(JV+13)
      GPHI=Q(JV+14)
      GPSI=Q(JV+15)
      SINPSI=SIN(GPSI*DEGRAD)
      COSPSI=COS(GPSI*DEGRAD)
      GU0=Q(JV+16)
      GV0=Q(JV+17)
      GSCU=Q(JV+18)
      GSCV=Q(JV+19)
      ICUT=Q(JV+22)
      NGVIEW=0
      IF (IVIEW.EQ.IDVIEW) GO TO 10
C
C             IVIEW is not current view bank
C
      MIGU=Q(JV+1)
      MITU=Q(JV+7)
      GO TO 20
C
C             IVIEW is current view bank
C
   10 MIGU=IGU
      MITU=ITU
C
C             Show graphic units
C
   20 IF (MIGU.EQ.0) GO TO 50
C
      DO 40 I=1,MIGU
C
C             Bit 1- 7 = LENGU
C             Bit 8-24 = Line attribute
C
         ISUM=Q(LQ(JV-1)+I)
*SG
         IFIL=IBITS(ISUM,13,3)
         IF(IFIL.EQ.0)THEN
            ILEN=IBITS(ISUM,0,7)
         ELSE
            ILEN=IBITS(ISUM,0,10)
         ENDIF
         CALL GWORK(2*ILEN)
         LINATT=IBITS(ISUM,0,24)
         IF(IFIL.EQ.0)THEN
            CALL MVBITS(0,0,7,LINATT,0)
         ELSE
            CALL MVBITS(0,0,10,LINATT,0)
         ENDIF
*SG
C
         IADD=Q(LQ(JV-2)+I)
C
         DO 30 J=1,ILEN
C
C             Zoom processing
C
            WS(J) = Q(LQ(JV-4)+IADD+J-1)*GZUA+GZUB+GZUC
            WS(ILEN+J) = Q(LQ(JV-5)+IADD+J-1)*GZVA+GZVB+GZVC
C
   30    CONTINUE
C
         MEMIDV=IDVIEW
         IF (IDVIEW.EQ.0) IDVIEW=-175
         CALL GDRAWV(WS(1),WS(ILEN+1),ILEN)
         IDVIEW=MEMIDV
         JV=LQ(JDRAW-IVIEW)
C
   40 CONTINUE
C
C             Show text units
C
   50 IF (MITU.EQ.0) GO TO 80
C
      DO 60 I=1,MITU
         IADD=Q(LQ(JV-3)+I)
         JVAD=LQ(JV-6)+IADD
         X=Q(JVAD)
         Y=Q(JVAD+1)
C
C             Zoom processing
C
         X=X*GZUA+GZUB+GZUC
         Y=Y*GZVA+GZVB+GZVC
C
         SIZE=Q(JVAD+2)*(GZUA+GZVA)*.5
         ANGLE=Q(JVAD+3)
C
C             Bit 1- 3 = Text line width (LWIDTH)
C             Bit 4-24 = Text attribute
C
         ISUM=Q(JVAD+4)
         LWIDTH=IBITS(ISUM,0,3)
         ITXATT=IBITS(ISUM,0,24)
C
         IOPT=Q(JVAD+5)
C
         NCHAR=Q(JVAD+6)
         CHTEXT=' '
         CALL UHTOC(IQ(JVAD+7),4,CHTEXT,NCHAR)
         MEMIDV=IDVIEW
         IF (IDVIEW.EQ.0) IDVIEW=-175
         CALL GDRAWT(X,Y,CHTEXT,SIZE,ANGLE,LWIDTH,IOPT)
         IDVIEW=MEMIDV
         JV=LQ(JDRAW-IVIEW)
C
   60 CONTINUE
C
      GO TO 80
C
   70 WRITE (CHMAIL,10000) IVIEW
      CALL GMAIL(0,0)
C
   80 CALL GDCOL(0)
      RETURN
10000 FORMAT(' *** GDSHOW ***: View ',I6,' does not exist')
      END
