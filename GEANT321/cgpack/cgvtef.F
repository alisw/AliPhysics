*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGVTEF(IFACE,XYZ,IVIS)
************************************************************************
*                                                                      *
*     Name: CGVTEF                                                     *
*     Author: E. Chernyaev                       Date:    16.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Test edge against face                                 *
*                                                                      *
*     References: CGHPLA, CGVTIN                                       *
*                                                                      *
*     Input: IFACE(*) - face                                           *
*            XYZ(3,*) - nodes coordinates                              *
*                                                                      *
*     Output: IVIS = -1 for unvisible edge                             *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgdelt.inc"
#include "geant321/cgcedg.inc"
*SG
      INTEGER   IFACE(*)
*SG
      REAL      XYZ(3,*),ABCD(4)
*-
      ERROR  = 2.*EESCR
      IVIS   = 0
      NEDGE  = IFACE(1)
      JF     = 2
      JT     = 0
      TMIN   = 99999.
      TMAX   =-99999.
*           C O M P U T E   A L L   I N T E R S E C T I O N   POINTS
      DO 140 NE=1,NEDGE
        N1     = IFACE(JF)
        N2     = IFACE(JF+1)
        JF     = JF + 2
        TY1    = ABC(1,2)*XYZ(1,N1) + ABC(2,2)*XYZ(2,N1) + ABC(3,2)
        TY2    = ABC(1,2)*XYZ(1,N2) + ABC(2,2)*XYZ(2,N2) + ABC(3,2)
        IF (ABS(TY1) .LE. ERROR) TY1 = 0.
        IF (ABS(TY2) .LE. ERROR) TY2 = 0.
        IF (TY1*TY2) 110,120,140
  110   TX1    = ABC(1,1)*XYZ(1,N1) + ABC(2,1)*XYZ(2,N1) + ABC(3,1)
        TX2    = ABC(1,1)*XYZ(1,N2) + ABC(2,1)*XYZ(2,N2) + ABC(3,1)
        T1     = (TX1*TY2-TX2*TY1) / (TY2-TY1)
        T2     = 0.
        GOTO 130
  120   IF (TY2 .NE. 0.) GOTO 121
        IF (TY1 .NE. 0.) GOTO 122
        TTT(1,JT+1) = ABC(1,1)*XYZ(1,N1) + ABC(2,1)*XYZ(2,N1) + ABC(3,1)
        TTT(2,JT+1) = 0
        TTT(1,JT+2) = ABC(1,1)*XYZ(1,N2) + ABC(2,1)*XYZ(2,N2) + ABC(3,1)
        TTT(2,JT+2) = 0
        IF (TTT(1,JT+1) .LT. TMIN)  TMIN = TTT(1,JT+1)
        IF (TTT(1,JT+1) .GT. TMAX)  TMAX = TTT(1,JT+1)
        IF (TTT(1,JT+2) .LT. TMIN)  TMIN = TTT(1,JT+1)
        IF (TTT(1,JT+2) .GT. TMAX)  TMAX = TTT(1,JT+1)
        GOTO 140
  121   T1     = ABC(1,1)*XYZ(1,N1) + ABC(2,1)*XYZ(2,N1) + ABC(3,1)
        IF (TY2 .LT. 0.) T2=-1.
        IF (TY2 .GT. 0.) T2= 1.
        GOTO 130
  122   T1     = ABC(1,1)*XYZ(1,N2) + ABC(2,1)*XYZ(2,N2) + ABC(3,1)
        IF (TY1 .LT. 0.) T2=-1
        IF (TY1 .GT. 0.) T2= 1
        GOTO 130
  130   JT     = JT + 1
        TTT(1,JT)= T1
        TTT(2,JT)= T2
        IF (TTT(1,JT) .LT. TMIN)  TMIN = TTT(1,JT)
        IF (TTT(1,JT) .GT. TMAX)  TMAX = TTT(1,JT)
  140   CONTINUE
*           C H E C K   I F   F A C E   O U T   O F   E D G E
      IF (JT .EQ. 0)                  GOTO 999
      IF (TMIN.GE.1.-TDEL .OR. TMAX.LE.0.+TDEL) GOTO 999
*           C O M P U T E   F A C E   P L A N E
*           C H E C K   I F   E D G E   B E F O R E  P L A N E
      CALL CGHPLA(IFACE,XYZ,ABCD)
      DA    = ABCD(1)*AA(1)+ABCD(2)*AA(2)+ABCD(3)*AA(3)+ABCD(4)
      DB    = ABCD(1)*BB(1)+ABCD(2)*BB(2)+ABCD(3)*BB(3)+ABCD(4)
      IF (DA.GE.-2.*EESCR .AND. DB.GE.-2.*EESCR)      GOTO 999
      IF (DA.LE.0. .AND. DB.LE.0.)              GOTO 150
      T0    = DA/(DA-DB)
      IF (DA.GE.0. .AND. TMAX.LE.T0+TDEL)       GOTO 999
      IF (DB.GE.0. .AND. TMIN.GE.T0-TDEL)       GOTO 999
*          S O R T   I N T E R S E C T I O N   P O I N T S
      IF (JT .EQ. 1)    GOTO 400
  150 DO 300 I=1,JT-1
        DO 200 J=I,1,-1
          IF (TTT(1,J+1) .GE. TTT(1,J))   GOTO 300
          T1         = TTT(1,J+1)
          T2         = TTT(2,J+1)
          TTT(1,J+1) = TTT(1,J)
          TTT(2,J+1) = TTT(2,J)
          TTT(1,J)   = T1
          TTT(2,J)   = T2
  200     CONTINUE
  300   CONTINUE
*
      DO 310 I=2,JT
        IF (TTT(2,I-1) .EQ. -2.)      GOTO 310
*
******* 02-FEB-1990  **** R.N.
        IF(.0001.LT.ABS(TTT(1,I)-TTT(1,I-1)))GOTO 310
******* 02-FEB-1990  **** R.N.
*
        IF (TTT(1,I) .NE. TTT(1,I-1)) GOTO 310
        IF (TTT(2,I) .EQ. TTT(2,I-1)) GOTO 305
        IF (TTT(2,I) .EQ. 0.)         GOTO 305
        IF (TTT(2,I-1) .EQ. 0.)       GOTO 305
        TTT(2,I)   = -2
        GOTO 310
  305   TTT(2,I-1) = -2.
        TTT(2,I)   = -2.
  310   CONTINUE
*
      J     = 0
      DO 320 I=1,JT
        IF (TTT(2,I) .EQ. -2.)        GOTO 320
        J     = J + 1
        TTT(1,J) = TTT(1,I)
  320   CONTINUE
      JT    = J
      IF (JT .EQ. 0)                  GOTO 999
*
  400 T1    = 0.
      T2    = 1.
      IF (DA .GT. 0.)           T1 = T0
      IF (DB .GT. 0.)           T2 = T0
      DO 500 I=1,JT,2
        IF (TTT(1,I+1) .LE. T1)   GOTO 500
        IF (TTT(1,I)   .GE. T2)   GOTO 999
        TMIN = TTT(1,I)
        TMAX = TTT(1,I+1)
        IF (TMIN .LT. T1)       TMIN = T1
        IF (TMAX .GT. T2)       TMAX = T2
        IF (TMAX-TMIN.LE.TDEL)  GOTO 500
        CALL CGVTIN(TMIN,TMAX,IVIS)
        IF (IVIS .LT. 0)        GOTO 999
  500   CONTINUE
  999 RETURN
      END
