*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGBTTT(WHAT,TMIN,TMAX,NT,NEDGE)
************************************************************************
*                                                                      *
*     Name: CGBTTT                                                     *
*     Author: E. Chernyaev                       Date:    14.03.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Prepare list of T-intervals                            *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: WHAT - flag: what is need ('GT','GE','LT','LE')           *
*            TMIN - min of T                                           *
*            TMAX - max of T                                           *
*              NT - number of T-points                                 *
*                                                                      *
*     Output: NEDGE - number of edges                                  *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgcedg.inc"
#include "geant321/cgdelt.inc"
      CHARACTER*2 WHAT
      CHARACTER*3 STATUS
*-
      TDEL   = 4.*EEWOR
      NE     = 0
      IF (NT .GE. 2)                    GOTO 100
*      IF (NT .EQ. 1)                    PRINT *,' CGBTTT: Number of T=1'
      IF (WHAT(1:1) .EQ. 'L')           GOTO 999
      TTT(1,1) = TMIN
      TTT(2,1) = TMAX
      ITTT(1)  = 0
      NE       = 1
      GOTO 999
*
**          S O R T    I N T E R S E C T I O N   P O I N T S
*
  100 DO 120 I=1,NT-1
        DO 110 J=I,1,-1
          IF (TTT(1,J+1) .GE. TTT(1,J)) GOTO 120
          T1         = TTT(1,J+1)
          T2         = TTT(2,J+1)
          TTT(1,J+1) = TTT(1,J)
          TTT(2,J+1) = TTT(2,J)
          TTT(1,J)   = T1
          TTT(2,J)   = T2
  110     CONTINUE
  120   CONTINUE
*
**         F I N D   S I G N   O F   1 - S T   P O I N T
*
  200 I1NEG  = 0
      I1POS  = 0
      I1ZBEG = 0
      I1ZEND = 0
      I1SIGN =-1
      T1     = TTT(1,1)
      DO 220 I=1,NT
        J      = I
        IF (TTT(1,I)-T1 .LE. TDEL)              GOTO 210
        IF (I1NEG .NE. I1POS)                   GOTO 240
        T1     = TTT(1,I)
        I1NEG  = 0
        I1POS  = 0
        I1ZBEG = 0
        I1ZEND = 0
  210   IF (TTT(2,I) .EQ. -1.)                  I1NEG  = 1
        IF (TTT(2,I) .EQ. +1.)                  I1POS  = 1
        IF (TTT(2,I) .EQ.  0.)                  I1ZBEG = 1
        IF (TTT(2,I) .EQ.  2.)                  I1ZEND = 1
  220   CONTINUE
      IF (I1NEG .EQ. I1POS)                     GOTO 250
  240 IF (I1NEG .EQ. 1)                         I1SIGN =-1
      IF (I1POS .EQ. 1)                         I1SIGN =+1
      IF (I1ZEND.EQ.1 .AND. I1ZBEG.EQ.0)        I1SIGN =-I1SIGN
*
  250 NNTT   = NT
      IF (WHAT(1:1) .EQ. 'L')   GOTO 300
      NNTT   = NT + 1
      TTT(1,NNTT) =+99999.
      TTT(2,NNTT) = I1SIGN
*
**          P R E P A R E   L I S T   O F   I N T E R V A L S
*
  300 T2     = -99999.
      I2NEG  = 0
      I2POS  = 0
      I2ZBEG = 0
      I2ZEND = 0
      IF (I1SIGN .LT. 0)        I2POS = 1
      IF (I1SIGN .GT. 0)        I2NEG = 1
      NE     = 0
      ITCUR  = 1
      STATUS = 'OUT'
*           S T A R T   O F   S E A R C H   F O R   N E X T   P O I N T
  310 T1     = T2
      I1NEG  = I2NEG
      I1POS  = I2POS
      I1ZBEG = I2ZBEG
      I1ZEND = I2ZEND
  320 IF (ITCUR .GT. NNTT)              GOTO 600
      T2     = TTT(1,ITCUR)
      I2NEG  = 0
      I2POS  = 0
      I2ZBEG = 0
      I2ZEND = 0
  330 IF (TTT(2,ITCUR) .EQ.-1.)         I2NEG  = 1
      IF (TTT(2,ITCUR) .EQ.+1.)         I2POS  = 1
      IF (TTT(2,ITCUR) .EQ. 0.)         I2ZBEG = 1
      IF (TTT(2,ITCUR) .EQ. 2.)         I2ZEND = 1
      ITCUR  = ITCUR + 1
      IF (ITCUR .GT. NNTT)              GOTO 400
      IF (TTT(1,ITCUR)-T2 .LE. TDEL)    GOTO 330
*           N E X T   P O I N T   I S   F O U N D E D
  400 T2     = (T2+TTT(1,ITCUR-1))/2.
      IF (I1ZBEG .NE. 0)                GOTO 410
      IF (I2ZBEG*I2ZEND .NE. 0)         GOTO 420
      IF (I2POS*I2NEG .NE. 0)           GOTO 530
      IF (I1SIGN.GT.0 .AND. I2POS.GT.0) GOTO 510
      IF (I1SIGN.LT.0 .AND. I2NEG.GT.0) GOTO 510
      IF (I1SIGN.GT.0 .AND. I2NEG.GT.0) GOTO 520
      IF (I1SIGN.LT.0 .AND. I2POS.GT.0) GOTO 520
      GOTO 310
*           B O U N D A R Y   E D G E
  410 NE        = NE + 1
      TTT(1,NE) = T1
      TTT(2,NE) = T2
      ITTT(NE)  = 1
      IF (I2ZEND .EQ. 0)                I2ZBEG = 1
  405 IF (I1SIGN.GT.0 .AND. I2POS.GT.0) STATUS = 'IN '
      IF (I1SIGN.LT.0 .AND. I2NEG.GT.0) STATUS = 'IN '
      IF (I1SIGN.GT.0 .AND. I2NEG.GT.0) STATUS = 'OUT'
      IF (I1SIGN.LT.0 .AND. I2POS.GT.0) STATUS = 'OUT'
      GOTO 310
*           V E R Y   S M A L L  B O U N D A R Y   E D G E
  420 IF (I2POS+I2NEG .EQ. 0)           GOTO 310
      I2ZBEG = 0
      I2ZEND = 0
      GOTO 405
*
**          C O M E   I N
*
  510 IF (STATUS .EQ. 'IN ')            GOTO 511
      IF (STATUS .EQ. 'OUT')            GOTO 512
  511 IF (WHAT(1:1) .EQ. 'L')           GOTO 550
      GOTO 310
  512 STATUS = 'IN '
      IF (WHAT(1:1) .EQ. 'G')           GOTO 550
      GOTO 310
*           C O M E   O U T
  520 IF (STATUS .EQ. 'IN ')            GOTO 521
      IF (STATUS .EQ. 'OUT')            GOTO 522
  521 STATUS = 'OUT'
      IF (WHAT(1:1) .EQ. 'L')           GOTO 550
      GOTO 310
  522 IF (WHAT(1:1) .EQ. 'L')           GOTO 550
      GOTO 310
*           C O M E   I N / O U T
  530 IF (STATUS .EQ. 'IN ')            GOTO 531
      IF (STATUS .EQ. 'OUT')            GOTO 532
  531 IF (WHAT(1:1) .EQ. 'L')           GOTO 550
      GOTO 310
  532 IF (WHAT(1:1) .EQ. 'G')           GOTO 550
      GOTO 310
*
  550 NE        = NE + 1
      TTT(1,NE) = T1
      TTT(2,NE) = T2
      ITTT(NE)  = 0
      GOTO 310
*
**          S K I P   B O U N D A R Y   E D G E S   (I F   N E E D)
*
  600 IF (NE .EQ. 0)            GOTO 999
      IF (WHAT(2:2) .EQ. 'E')   GOTO 700
      K      = 0
      DO 610 I=1,NE
        IF (ITTT(I) .EQ. 1)     GOTO 610
        K      = K + 1
        TTT(1,K) = TTT(1,I)
        TTT(2,K) = TTT(2,I)
        ITTT(K)  = ITTT(I)
  610   CONTINUE
      NE     = K
*
  700 IF (NE .EQ. 0)            GOTO 999
      IF (TMIN .GT. TMAX)       GOTO 999
      K      = 0
      DO 710 I=1,NE
        IF (TTT(2,I) .LE. TMIN) GOTO 710
        IF (TTT(1,I) .GE. TMAX) GOTO 720
        T1   = TTT(1,I)
        T2   = TTT(2,I)
        IF (T1 .LT. TMIN)       T1 = TMIN
        IF (T2 .GT. TMAX)       T2 = TMAX
        IF (T2-T1 .LT. TDEL)    GOTO 710
        K      = K + 1
        TTT(1,K) = T1
        TTT(2,K) = T2
        ITTT(K)  = ITTT(I)
  710   CONTINUE
  720 NE     = K
*
  999 NEDGE  = NE
      RETURN
      END
