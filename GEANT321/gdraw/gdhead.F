*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:23  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.26  by  S.Giani
*-- Author :
      SUBROUTINE GDHEAD(ISEL,NAME,CHRSIZ)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       This routine draws a graphics header.                    *
C.    *       NAME is an hollerith string terminated by a dollar sign. *
C.    *                                                                *
C.    *       Different options can be chosen by the parameter ISEL :  *
C.    *       =0 to have only the header lines                         *
C.    *       =xxxxx1 to add the text NAME centered on top of header   *
C.    *       =xxxx1x to add global detector name (first vol.) on left *
C.    *       =xxx1xx to add date on right                             *
C.    *       =xx1xxx to select thick characters for NAME (i.e. with   *
C.    *               larger line width)                               *
C.    *       =x1xxxx to add the text 'Event NR x' on top of header    *
C.    *       =1xxxxx to add the text 'Run NR x' on top of header      *
C.    *                                                                *
C.    *       A negative value of ISEL means special headers, i.e.     *
C.    *       directly used by GEANT3 routines :                       *
C.    *       =-1 : GDSPEC header                                      *
C.    *       =-3 : VIEWXY header                                      *
C.    *       =-4 : VIEWYZ header                                      *
C.    *       =-5 : VIEWXZ header                                      *
C.    *                                                                *
C.    *    ==>Called by : <USER>, <GXINT>, GDSPEC, GDTREE              *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcflag.inc"
#include "geant321/gcdraw.inc"
#include "geant321/gctime.inc"
C
      DIMENSION UH(5),VH(5),UHH(2),VHH(2)
      DIMENSION UH1(5),VH1(5),UHH1(2),VHH1(2)
      CHARACTER*9 CHDATE
      CHARACTER*4 CHTEXT
      CHARACTER*(*) NAME
C
      SAVE UH1,VH1,UHH1,VHH1
C
      DATA UH1/0.1,0.1,19.9,19.9,0.1/
      DATA VH1/19.9,0.1,0.1,19.9,19.9/
      DATA UHH1/0.1,19.9/
      DATA VHH1/18.5,18.5/
C.
C.    ------------------------------------------------------------------
C.
C             Normalize to PLTRNX,PLTRNY
C
      DO 10 I=1,5
         UH(I)=UH1(I)*PLTRNX/20.
         VH(I)=VH1(I)*PLTRNY/20.
   10 CONTINUE
      DO 20 I=1,2
         UHH(I)=UHH1(I)*PLTRNX/20.
         VHH(I)=VHH1(I)*PLTRNY/20.
   20 CONTINUE
C
      IF (ISEL.GE.0) THEN
         ISEL5=ISEL/100000
         ISEL5=ISEL5-ISEL5/2*2
         ISEL4=ISEL/10000
         ISEL4=ISEL4-ISEL4/2*2
         ISEL3=ISEL/1000
         ISEL3=ISEL3-ISEL3/2*2
         ISEL2=ISEL/100
         ISEL2=ISEL2-ISEL2/2*2
         ISEL1=ISEL/10
         ISEL1=ISEL1-ISEL1/2*2
         ISEL0=ISEL/1
         ISEL0=ISEL0-ISEL0/2*2
      ENDIF
C
      CALL GDCOL(1)
      CALL GDLW(2)
      CALL GDRAWV(UH,VH,5)
      IF (ISEL.EQ.0) THEN
         CALL GDLW(0)
         GO TO 30
      ENDIF
      CALL GDRAWV(UHH,VHH,2)
      CALL GDLW(0)
      IF (ISEL.LT.0.OR.ISEL1.NE.0) THEN
         CALL UHTOC(IQ(JVOLUM+1),4,CHTEXT,4)
         XTEXT=1.*PLTRNX/20.
         YTEXT=19.*PLTRNY/20.
         CSIZE=0.3*MIN(PLTRNX,PLTRNY)/20.
         CALL GDRAWT(XTEXT,YTEXT,CHTEXT,CSIZE,0.,1,-1)
      ENDIF
      IF (ISEL3.NE.0) THEN
         LWDTH=2
      ELSE
         LWDTH=1
      ENDIF
      IF (ISEL0.NE.0) THEN
         CALL GDCOL(4)
         XTEXT=10.*PLTRNX/20.
         YTEXT=19.*PLTRNY/20.
         CSIZE=CHRSIZ*MIN(PLTRNX,PLTRNY)/20.
         CALL GDRAWT(XTEXT,YTEXT,NAME,CSIZE,0.,LWDTH,0)
         CALL GDCOL(1)
      ENDIF
      IF (ISEL.GT.0) THEN
         IF (ISEL4.NE.0) THEN
            XTEXT=7.*PLTRNX/20.
            YTEXT=18.7*PLTRNY/20.
            CSIZE=0.3*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,'EVENT',CSIZE,0.,LWDTH,-1)
            XTEXT=9.*PLTRNX/20.
            CALL GDRAWT(XTEXT,YTEXT,'NR',CSIZE,0.,LWDTH,-1)
*
            CHDATE=' '
            WRITE (CHDATE,'(I6)') IDEVT
            XTEXT=12.*PLTRNX/20.
            YTEXT=18.7*PLTRNY/20.
            CSIZE=0.3*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,CHDATE,CSIZE,0.,LWDTH,1)
         ENDIF
         IF (ISEL5.NE.0) THEN
            XTEXT=7.*PLTRNX/20.
            YTEXT=19.2*PLTRNY/20.
            CSIZE=0.3*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,'RUN',CSIZE,0.,LWDTH,-1)
            XTEXT=9.*PLTRNX/20.
            CALL GDRAWT(XTEXT,YTEXT,'NR',CSIZE,0.,LWDTH,-1)
*
            CHDATE=' '
            WRITE (CHDATE,'(I6)') IDRUN
            XTEXT=12.*PLTRNX/20.
            YTEXT=19.2*PLTRNY/20.
            CSIZE=0.3*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,CHDATE,CSIZE,0.,LWDTH,1)
         ENDIF
      ELSE
         IF (ISEL.EQ.-1) THEN
            CALL GDCOL(4)
            XTEXT=8.5*PLTRNX/20.
            YTEXT=19.*PLTRNY/20.
            CSIZE=0.6*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,NAME,CSIZE,0.,2,1)
            CALL GDCOL(1)
            XTEXT=9.*PLTRNX/20.
            YTEXT=19.*PLTRNY/20.
            CSIZE=0.4*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT,'<SPECIFICATIONS',CSIZE,0.,1,-1)
         ELSE IF (ISEL.EQ.-3) THEN
            CALL GDCOL(4)
            XTEXT=10.*PLTRNX/20.
            YTEXT=19.*PLTRNY/20.
            CSIZE=0.5*MIN(PLTRNX,PLTRNY)/20.
            CALL GDRAWT(XTEXT,YTEXT, 'C<ROSS - >S<ECTION> (<Z >=< >0)',
     +      CSIZE,0.,2,0)
            CALL GDCOL(1)
         ELSE IF (ISEL.EQ.-4) THEN
            CALL GDCOL(4)
            XTEXT=10.*PLTRNX/20.
            YTEXT=19.*PLTRNY/20.
            CSIZE=0.5*MIN(PLTRNX,PLTRNY)/20.
          CALL GDRAWT(XTEXT,YTEXT,
     +      'V<ERTICAL> S<ECTION> (<X >=< >0.5)',CSIZE,0.,2,0)
            CALL GDCOL(1)
         ELSE IF (ISEL.EQ.-5) THEN
            CALL GDCOL(4)
            XTEXT=10.*PLTRNX/20.
            YTEXT=19.*PLTRNY/20.
            CSIZE=0.5*MIN(PLTRNX,PLTRNY)/20.
          CALL GDRAWT(XTEXT,YTEXT,
     +      'H<ORIZONTAL> S<ECTION> (<Y >=< >5)',CSIZE,0.,2,0)
            CALL GDCOL(1)
         ENDIF
      ENDIF
      IF (ISEL.LT.0.OR.ISEL2.NE.0) THEN
         IYEAR  = (IGDATE)/10000
         IMONTH = (IGDATE-IYEAR*10000)/100
         IDAY   = (IGDATE-IYEAR*10000-IMONTH*100)
         CHDATE=' '
         WRITE (CHDATE,'(I2,''/'',I2,''/'',I2)') IDAY,IMONTH,IYEAR
         XTEXT=19.*PLTRNX/20.
         YTEXT=19.*PLTRNY/20.
         CSIZE=0.25*MIN(PLTRNX,PLTRNY)/20.
         CALL GDRAWT(XTEXT,YTEXT,CHDATE,CSIZE,0.,1,1)
      ENDIF
   30 CALL GDCOL(0)
      END
