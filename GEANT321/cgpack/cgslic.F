*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGSLIC(A,SLIC,NMAX,B)
************************************************************************
*                                                                      *
*     Name: CGSLIC                                                     *
*     Author: E. Chernyaev                       Date:    20.04.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Make slice of CG-object                                *
*                                                                      *
*     References: CGMNMX                                               *
*                                                                      *
*     Input:  A(*) - CG-object                                         *
*          SLIC(4) - slicing plane                                     *
*            NMAX  - max size of B-array                               *
*                                                                      *
*     Output:  B - resulting CG-object                                 *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      PARAMETER (NWORK=LCGHEA+6*(LCGFAC+4*LCGEDG))
      REAL      A(*),B(*),SLIC(4),ABCD(4)
      REAL      RMN(3),RMX(3),RMID(3),W(NWORK),XYZ(3,8)
*-
      B(1)   = 0.
*           T E S T   P A R A M E T E R S    C O R R E C T N E S S
      IF (NMAX .LT. LCGHEA)     GOTO 999
      S      = SQRT(SLIC(1)*SLIC(1)+SLIC(2)*SLIC(2)+SLIC(3)*SLIC(3))
      IF (S .LE. EEWOR)         GOTO 999
      ABCD(1)= SLIC(1) / S
      ABCD(2)= SLIC(2) / S
      ABCD(3)= SLIC(3) / S
      ABCD(4)= SLIC(4) / S
      CALL CGSCOP(1,A,RMN,RMX)
*      CALL CGMNMX(A,RMN,RMX)
      IF (RMN(1) .GT. RMX(1))   GOTO 999
*           M I N - M A X    T E S T
      NFACE  = A(KCGNF)
      IF (NFACE .EQ. 0)         GOTO 998
      NPOS   = 0
      NNEG   = 0
      J      = LCGHEA
      DO 120 NF=1,NFACE
        NEDGE  = A(J+KCGNE)
        J      = J + LCGFAC
        DO 110 NE=1,NEDGE
          DIST = ABCD(1)*A(J+KCGX1)+ABCD(2)*A(J+KCGY1)
     +          +ABCD(3)*A(J+KCGZ1) + ABCD(4)
          IF (DIST .GT.-EEWOR)    NPOS = NPOS + 1
          IF (DIST .LT.+EEWOR)    NNEG = NNEG + 1
          J      = J + LCGEDG
  110   CONTINUE
  120 CONTINUE
      IF (NPOS .EQ. 0)            GOTO 998
      IF (NNEG .EQ. 0)            CALL CGCOPY(A,NMAX,B)
      IF (NNEG .EQ. 0)            GOTO 999
*           P R E P A R E   S L I C I N G   O B J E C T
      K      = 1
      IF (ABS(ABCD(2)) .GT. ABS(ABCD(1)))       K = 2
      IF (ABS(ABCD(3)) .GT. ABS(ABCD(K)))       K = 3
      RMID(1) = (RMN(1)+RMX(1)) / 2.
      RMID(2) = (RMN(2)+RMX(2)) / 2.
      RMID(3) = (RMN(3)+RMX(3)) / 2.
      RX      = RMX(1) - RMN(1)
      RY      = RMX(2) - RMN(2)
      RZ      = RMX(3) - RMN(3)
*
      GOTO (210,220,230),K
  210 XYZ(2,1) = RMID(2) + RY
      XYZ(3,1) = RMID(3) + RZ
      XYZ(2,2) = RMID(2) - RY
      XYZ(3,2) = RMID(3) + RZ
      XYZ(2,3) = RMID(2) - RY
      XYZ(3,3) = RMID(3) - RZ
      XYZ(2,4) = RMID(2) + RY
      XYZ(3,4) = RMID(3) - RZ
      DO 215 I=1,4
        XYZ(1,I) = -(ABCD(2)*XYZ(2,I)+ABCD(3)*XYZ(3,I)+ABCD(4))/ABCD(1)
  215   CONTINUE
      GOTO 250
*
  220 XYZ(1,1) = RMID(1) + RX
      XYZ(3,1) = RMID(3) + RZ
      XYZ(1,2) = RMID(1) - RX
      XYZ(3,2) = RMID(3) + RZ
      XYZ(1,3) = RMID(1) - RX
      XYZ(3,3) = RMID(3) - RZ
      XYZ(1,4) = RMID(1) + RX
      XYZ(3,4) = RMID(3) - RZ
      DO 225 I=1,4
        XYZ(2,I) = -(ABCD(1)*XYZ(1,I)+ABCD(3)*XYZ(3,I)+ABCD(4))/ABCD(2)
  225   CONTINUE
      GOTO 250
*
  230 XYZ(1,1) = RMID(1) + RX
      XYZ(2,1) = RMID(2) + RY
      XYZ(1,2) = RMID(1) - RX
      XYZ(2,2) = RMID(2) + RY
      XYZ(1,3) = RMID(1) - RX
      XYZ(2,3) = RMID(2) - RY
      XYZ(1,4) = RMID(1) + RX
      XYZ(2,4) = RMID(2) - RY
      DO 235 I=1,4
        XYZ(3,I) = -(ABCD(1)*XYZ(1,I)+ABCD(2)*XYZ(2,I)+ABCD(4))/ABCD(3)
  235   CONTINUE
      GOTO 250
*
  250 IF (ABCD(K) .GT. 0)     S = -3.
      IF (ABCD(K) .LT. 0)     S = +3.
      RX = 0.
      RY = 0.
      RZ = 0.
      IF (K .EQ. 1)           RX = S * (RMX(1)-RMN(1))
      IF (K .EQ. 2)           RY = S * (RMX(2)-RMN(2))
      IF (K .EQ. 3)           RZ = S * (RMX(3)-RMN(3))
      DO 255 I=1,4
        XYZ(1,I+4) = XYZ(1,I) + RX
        XYZ(2,I+4) = XYZ(2,I) + RY
        XYZ(3,I+4) = XYZ(3,I) + RZ
  255   CONTINUE
      CALL CGBOX(XYZ,4,4,NWORK,W)
      CALL CGCEV(1,W)
*           M A K E   S L I C E
      CALL CGSUB(A,W,NMAX,B)
      GOTO 999
*
  998 B(KCGSIZ) = LCGHEA
      B(KCGATT) = 0.
      B(KCGNF)  = 0.
  999 RETURN
      END
