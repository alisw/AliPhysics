*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.28  by  S.Giani
*-- Author :
      SUBROUTINE GFUPAR(JVO,JIN,NPAR,NATT,PAR,ATT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *    ROUTINE TO FETCH USER PARAMETERS AND ATTRIBUTES FOR         *
C.    *    THE VOLUME FROM VOLUME ADDRESS JVO AND IF GSPOSP TYPE       *
C.    *    VERSION DEFINED BY IN ADDRESS JIN.                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>                                       *
C.    *         Author  P.Zanarini  *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcunit.inc"
#include "geant321/gconsp.inc"
      DIMENSION PAR(*),ATT(*)
C.
C.               ---------------------------------------
C.
      NPAR=Q(JVO+5)
      NATT=Q(JVO+6)
      ISH=Q(JVO+2)
      JPAR=JVO+7
      JATT=JVO+7+NPAR
C
      IF(NPAR.GT.0) GO TO 10
      IF(NPAR.NE.0) GO TO 900
      IF(JIN.LE.0) GO TO 910
C
      JPAR=JIN+9
      NPAR=Q(JPAR)
      JPAR=JPAR+1
C
   10 CONTINUE
C
C             NPAR : 30 ---> 12
C
      IF(ISH.EQ.28) NPAR=12
C
      IF(NPAR.LE.0) GO TO 920
      CALL UCOPY(Q(JPAR),PAR,NPAR)
      CALL UCOPY(Q(JATT),ATT,NATT)
C
      IF(ISH.NE.4) GO TO 20
C
C              TRAPEZOID
C
      PH=90.0
      IF(PAR(2).NE.0.)PH=ATAN2(PAR(3),PAR(2))*RADDEG
      TT=PAR(2)**2+PAR(3)**2
      IF(TT.GT.0.0) TT=SQRT(TT)
      PAR(2)=ATAN(TT)*RADDEG
      IF(PH.LT.0.0) PH=PH+360.0
      PAR(3)=PH
      PAR(7)=ATAN(PAR(7))*RADDEG
      IF(PAR(7).GT.90.0) PAR(7)=PAR(7)-180.0
      PAR(11)=ATAN(PAR(11))*RADDEG
      IF(PAR(11).GT.90.0) PAR(11)=PAR(11)-180.0
C
      GO TO 99
   20 CONTINUE
C
      IF(ISH.NE.10) GO TO 99
C
C              PARALLELEPIPED.
C
      PH=90.0
      IF(PAR(5).NE.0.)PH=ATAN2(PAR(6),PAR(5))*RADDEG
      TT=PAR(5)**2+PAR(6)**2
      IF(TT.GT.0.0) TT=SQRT(TT)
      PAR(4)=ATAN(PAR(4))*RADDEG
      IF(PAR(4).GT.90.0) PAR(4)=PAR(4)-180.0
      PAR(5)=ATAN(TT)*RADDEG
      IF(PH.LT.0.0) PH=PH+360.0
      PAR(6)=PH
C
      GO TO 99
C
  900 CONTINUE
C
      WRITE(CHMAIL,1000) NPAR
      CALL GMAIL(0,0)
 1000 FORMAT(' **** ERROR IN GFUPAR - NPAR FROM JVO =',I5,' LE 0 ****')
C
      GO TO 99
  910 CONTINUE
C
      WRITE(CHMAIL,1010) JIN
      CALL GMAIL(0,0)
 1010 FORMAT(' **** ERROR IN GFUPAR - NPAR FROM VOL 0 AND JIN ='
     +,I12,' IS LE 0 ****')
C
      GO TO 99
  920 CONTINUE
C
      WRITE(CHMAIL,1020) NPAR
      CALL GMAIL(0,0)
 1020 FORMAT(' **** ERROR IN GFUPAR - NPAR FROM JIN =',I5,' LE 0 ****')
C
   99 CONTINUE
      END
