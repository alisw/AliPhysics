*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE CONT(F,NUMBCS,CONTUR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PURPOSE                                                            C
C     TO PRODUCE A CONTOUR PLOT BY MEANS OF THE PRINTER                C
C   USAGE                                                              C
C     THIS DEPENDS ON WHETHER THE FUNCTION TO BE PLOTTED HAS REAL      C
C     ARGUMENTS OR INTEGER ARGUMENTS.                                  C
C     THE FUNCTION HAS REAL ARGUMENTS, SUPPOSE THE FUNCTION HAS        C
C     THE NAME Z.                                                      C
C     EXTERNAL Z                                                       C
C     DIMENSION CONTUR(10)                                             C
C     CALL PAPER....                                                   C
C     CALL NAMES....                                                   C
C     CALL FRAME....                                                   C
C     CONTUR(1)=....                                                   C
C     CONTUR(2)=....                                                   C
C     .                                                                C
C     .                                                                C
C     .                                                                C
C     CALL CONT(Z,NUMBCS,CONTUR)                                       C
C     IF YOU DO NOT KNOW WHAT SORT OF VALUES TO ASSIGN TO THE ARRAY    C
C     CONTUR, YOU CAN CALL THE SUBROUTINE FINDEM WHICH WILL HELP YOU.  C
C     CALL FINDEM(Z,NUMBCS,CONTUR)                                     C
C                                                                      C
C     F HAS INTEGER ARGUMENTS---                                       C
C     EXTERNAL FREARG                                                  C
C     DIMENSION CONTUR(10)                                             C
C     CALL PAPER....                                                   C
C     CALL NAMES....                                                   C
C     CALL FRAME....                                                   C
C     CONTUR(1)=....                                                   C
C     CONTUR(2)=....                                                   C
C     .                                                                C
C     .                                                                C
C     .                                                                C
C     CALL SETUP....                                                   C
C     CALL CONT(FREARG,NUMBCS,CONTUR)                                  C
C     AGAIN, YOU MAY USE FINDEM TO ASSIGN VALUES TO THE ARRAY CONTUR.  C
C     CALL FINDEM(FREARG,NUMBCS,CONTUR)                                C
C     NOTE THAT THE FIRST PARAMETER OF CONT AND FINDEM MUST HAVE       C
C     THE NAME FREARG.                                                 C
C     THE FUNCTION TO BE PLOTTED MUST HAVE THE NAME FINARG.            C
C   PARAMETERS                                                         C
C     F      - THE FUNCTION TO BE PLOTTED. F MUST BE A REAL FUNCTION   C
C              WITH TWO REAL PARAMETERS. IF THE FUNCTION YOU WANT TO   C
C              PLOT HAS INTEGER ARGUMENTS, THIS FUNCTION SHOULD BE     C
C              NAMED FINARG AND THE FIRST ACTUAL PARAMETER IN THE      C
C              CALL TO CONT SHOULD BE FREARG.                          C
C              (FREARG IS A FUNCTION WHICH CALLS FINARG FOUR TIMES     C
C              AND PERFORMS A TWO-WAY LINEAR INTERPOLATION ON THESE    C
C              FOUR VALUES.)                                           C
C     NUMBCS - THE NUMBER OF CONTOURS TO BE PLOTTED. NUMBS MUST BE     C
C              AT LEAST 1 AND AT MOST 10.                              C
C     CONTUR - AN ARRAY CONTAINING THE FUNCTION VALUES ON THE VARIOUS  C
C              CONTOURS. THESE FUNCTION VALUES DO NOT HAVE TO BE       C
C              EQUIDISTANT AND NOT EVEN MONOTONOUS.                    C
C              THE NUMBER 0 IS PRINTED WHERE F=CONTUR(1),              C
C              THE NUMBER 1 IS PRINTED WHERE F=CONTUR(2),              C
C              ETC.                                                    C
C   REMARKS                                                            C
C     THE CHOSEN NAME FOR F MUST APPEAR IN AN EXTERNAL STATEMENT IN    C
C     THE PROGRAM CALLING CONT.                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      EXTERNAL F
      DIMENSION CONTUR(10),SAVEC(10),SAVES(10)
      DIMENSION XVALUS(12)
      COMMON /J509C1/XBL,YBL,XUR,YUR,NX(6),NY(6),LETT1,LETT2,LINE1,LINE2
      REAL FU(133),FB(133)
      INTEGER LETTRI,MINUS,STAR,BLANK,SYMBOL(10),LINEIM(132)
      DATA LETTRI,MINUS,STAR,BLANK/1HI,1H-,1H*,1H /
#if !defined(CERNLIB_F90)
      DATA SYMBOL(1),SYMBOL(2),SYMBOL(3),SYMBOL(4)/1H0,1H1,1H2,1H3/
      DATA SYMBOL(5),SYMBOL(6),SYMBOL(7),SYMBOL(8)/1H4,1H5,1H6,1H7/
      DATA SYMBOL(9),SYMBOL(10)/1H8,1H9/
      DATA LSLASH,LDOT/1H/ , 1H. /
C=======================================================================
#endif
#if defined(CERNLIB_F90)
      INTEGER, DIMENSION(10) :: SYMBOL = (/ transfer('0   ', 0),        &
      transfer('1   ', 0), transfer('2   ', 0), transfer('3   ', 0),    &
      transfer('4   ', 0), transfer('5   ', 0), transfer('6   ', 0),    &
      transfer('7   ', 0), transfer('8   ', 0), transfer('9   ', 0) /)
      INTEGER :: LSLASH = transfer('/   ', 0), LDOT = transfer('.   ',0)
!=======================================================================
#endif
      IF (1.LE.NUMBCS.AND.NUMBCS.LE.10) GOTO 10
      WRITE (6, 1010) NUMBCS
      IF(1.LE.NUMBCS .AND. NUMBCS.LE.10)GO TO 10
      WRITE(6,1010)NUMBCS
      RETURN
   10 CONTINUE
C
C         INITIALIZE COMMON BLOCK IF NOT ALREADY DONE
      CALL J509BD
      LETT0 = LETT1 - 1
      LINE0 = LINE1 - 1
      LETT3 = LETT2 + 1
      LINE3 = LINE2 + 1
C         SAVE CONTOUR VALUES AND ORDER THEM ASCENDING
      CALL UCOPY(CONTUR,SAVEC,NUMBCS)
      CALL UCOPY(SYMBOL,SAVES,NUMBCS)
      CALL ORDRE2(CONTUR,SYMBOL,NUMBCS)
C
C   THE ARRAY LINEIM IS A BUFFER FOR A LINE IMAGE.                     C
C         PRINT LINE ALONG TOP OF FRAME
      WRITE(6,1030)NY
      LINEIM(LETT0)=STAR
      DO 30 LETT=LETT1,LETT2
   30 LINEIM(LETT)=MINUS
      LINEIM(LETT3)=STAR
      WRITE(6,1060)YUR,(LINEIM(LETT),LETT=LETT0,LETT3)
      MIDDLE=(LINE1+LINE2)/2
      DLINE = LINE2-LINE1+1
      DLETT = LETT2-LETT1+1
      DXVAL = (XUR-XBL) / DLETT
      DYVAL = (YUR-YBL) / DLINE
      YVAL = YUR
      XVAL = XBL
      DO 40 LETT=LETT1,LETT3
      FU(LETT) = F(XVAL,YVAL)
   40 XVAL = XVAL + DXVAL
CCC
      MAXCON = 1
      MINCON = 10
      YVAL = YUR
CCC       LOOP OVER LINES STARTS HERE
C
      DO 90 LINE=LINE1,LINE2
      YVAL = YVAL - DYVAL
      XVAL = XBL-DXVAL
      DO 50 LETT=LETT1,LETT3
      XVAL = XVAL + DXVAL
   50 FB(LETT) = F(XVAL,YVAL)
      DO 70 LETT=LETT1,LETT2
      FUL=FU(LETT)
      FUR=FU(LETT+1)
      FBL=FB(LETT)
      FBR=FB(LETT+1)
C   FUL IS THE FUNCTION VALUE AT THE POINT (LETT-1/2,LINE-1/2)         C
C   FUR IS THE FUNCTION VALUE AT THE POINT (LETT+1/2,LINE-1/2)         C
C   FBL IS THE FUNCTION VALUE AT THE POINT (LETT-1/2,LINE+1/2)         C
C   FBR IS THE FUNCTION VALUE AT THE POINT (LETT+1/2,LINE+1/2)         C
C                                                                      C
C   FIND OUT IF THERE IS A CONTOUR VALUE WHICH FALLS BETWEEN THE       C
C   BIGGEST AND THE SMALLEST OF THE FUNCTION VALUES AT THE FOUR        C
C   CORNERS OF THE PRESENT PRINT POSITION (LETT,LINE).                 C
C   IF SO, PRINT THE CORRESPONDING SYMBOL.                             C
C   IF NOT, PRINT A BLANK AT THE POINT (LETT,LINE).                    C
      BIG  =MAX(FUL,FUR,FBL,FBR)
      SMALL=MIN(FUL,FUR,FBL,FBR)
      LINEIM(LETT)=BLANK
      DO 60 NUMBC=1,NUMBCS
      IF (BIG   .LT. CONTUR(NUMBC))  GO TO 70
      IF (SMALL .GT. CONTUR(NUMBC))  GO TO 60
      LINEIM(LETT) = SYMBOL(NUMBC)
      IF (NUMBC .GT. MAXCON)  MAXCON = NUMBC
      IF (NUMBC .LT. MINCON)  MINCON = NUMBC
   60 CONTINUE
   70 CONTINUE
      IF (MOD(LINE,10) .EQ. LINE0)  GO TO 73
C         PRINT LINE (NORMAL LINE)
      LINEIM(LETT0) = LETTRI
      LINEIM(LETT3) = LETTRI
      WRITE(6,1050)      (LINEIM(LETT),LETT=LETT0,LETT3)
      GO TO 76
C         PRINT LINE (EVEN TENTH LINE)
   73 LINEIM(LETT0) = LDOT
      LINEIM(LETT3) = LDOT
      WRITE(6,1040) YVAL,(LINEIM(LETT),LETT=LETT0,LETT3)
   76 CONTINUE
      DO 80 LETT=LETT1,LETT3
   80 FU(LETT)=FB(LETT)
   90 CONTINUE
C         PRINT LINE ALONG BOTTOM OF FRAME
      LINEIM(LETT0)=STAR
      DO 100 LETT=LETT1,LETT2
  100 LINEIM(LETT)=MINUS
      LINEIM(LETT3)=STAR
      WRITE(6,1060)YBL,(LINEIM(LETT),LETT=LETT0,LETT3)
C         PRINT LABELS ON X-AXIS EVERY 10 COLIMNS
      LETT2P = LETT2 + 10
      DO 119 LETT= 2, LETT2P
  119 LINEIM(LETT) = BLANK
      X10COL = 10.*DXVAL
      ITENMX = 11
      LETT = LETT0
      DO 129 ITEN= 1, ITENMX
      XVALUS(ITEN) = XBL + X10COL*(ITEN-1)
      LINEIM(LETT) = LSLASH
      LETT = LETT + 10
      IF (XVALUS(ITEN) .GT. (XUR-5.*DXVAL)) GO TO 139
  129 CONTINUE
      ITEN = ITENMX
  139 CONTINUE
      WRITE(6,1090)  (LINEIM(LETT),LETT=2,LETT2P)
      WRITE(6,1100) (XVALUS(II),II=1,ITEN)
C         PRINT VARIABLE NAME FOR X-AXIS
      WRITE(6,1110)NX,DXVAL
      IF (MAXCON .LT. MINCON)  GO TO 125
      WRITE(6,1120) (SYMBOL(NN),CONTUR(NN),NN=MINCON,MAXCON)
      GO TO 135
  125 CONTINUE
      WRITE(6,1130)
  135 CONTINUE
C         RESTORE INITIAL CONTOUR ORDERING
      CALL UCOPY(SAVES,SYMBOL,NUMBCS)
      CALL UCOPY(SAVEC,CONTUR,NUMBCS)
      RETURN
 1010 FORMAT(' HOW MANY CONTOURS DID YOU WANT---',I10,'---(CONT)')
 1030 FORMAT('1',120A1)
 1040 FORMAT(1X,F8.3,1X,122A1)
 1050 FORMAT(10X,122A1)
 1060 FORMAT(1X,F8.3,'-',122A1)
 1090 FORMAT(1X,131A1)
 1100 FORMAT(4X,12(1X,F9.3))
 1110 FORMAT(20X,6A1,20X,'ONE COLUMN =',E10.3,10X,
     *    2('CONTOUR FCN VALUE     '))
 1120 FORMAT(80X,A1,F14.5,7X,A1,F14.5)
 1130 FORMAT(80X,'NO CONTOURS FOUND'/)
      END
