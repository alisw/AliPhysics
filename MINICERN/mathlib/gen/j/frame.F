*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:50  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      SUBROUTINE FRAME(XBLCOR,YBLCOR,XURCOR,YURCOR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PURPOSE                                                            C
C     TO SET UP THE PLOT FRAME FOR THE CONTOUR PLOTTING PROGRAM CONT.  C
C   USAGE                                                              C
C     CALL FRAME(XBLCOR,YBLCOR,XURCOR,YURCOR)                          C
C     THIS STATEMENT SHOULD BE EXECUTED BEFORE CONT IS CALLED.         C
C   PARAMETERS                                                         C
C     XBLCOR - THE X-COORDINATE OF THE BOTTOM LEFT CORNER OF THE       C
C              PLOT FRAME                                              C
C     YBLCOR - THE Y-COORDINATE OF THE BOTTOM LEFT CORNER OF THE       C
C              PLOT FRAME                                              C
C     XURCOR - THE X-COORDINATE OF THE UPPER RIGHT CORNER OF THE       C
C              PLOT FRAME                                              C
C     YURCOR - THE Y-COORDINATE OF THE UPPER RIGHT CORNER OF THE       C
C              PLOT FRAME                                              C
C
      COMMON /J509C1/XBL,YBL,XUR,YUR,NX(6),NY(6),LETT1,LETT2,LINE1,LINE2
      CALL J509BD
C
      XBL = XBLCOR
      YBL=YBLCOR
      XUR=XURCOR
      YUR=YURCOR
      IF(XUR.LE.XBL.OR.YUR.LE.YBL) GOTO 1
      RETURN
    1 WRITE(6,100)XBL,YBL,XUR,YUR
  100 FORMAT(' ERROR IN FRAME - EXIT CALLED'/' XBLCOR = ',F8.3,5X
     +  ,'YBLCOR = ',F8.3,5X,'XURCOR = ',F8.3,5X,'YURCOR = ',F8.3)
      STOP
      END
