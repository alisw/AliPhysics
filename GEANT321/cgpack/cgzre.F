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
      SUBROUTINE CGZRE(K,CG,LCG)
************************************************************************
*                                                                      *
*     Name: CGZROT                                                     *
*     Author: E. Chernyaev                       Date:    01.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Rotate edge about Z-axis                               *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: K  - number for step                                      *
*            CG - CG-object                                            *
*                                                                      *
*     Output: LCG - last index in CG-object                            *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cgcaan.inc"
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION  X1,Y1,H1,X2,Y2,H2,AW,AA,BB,CC,DD,S
      DOUBLE PRECISION  D1,D2,DX1,DY1,DH1,DX2,DY2,DH2
#endif
      REAL      CG(*)
*-
      J      = 0
      D1     = ABS(XYHA(1,1)) + ABS(XYHA(2,1))
      D2     = ABS(XYHA(1,2)) + ABS(XYHA(2,2))
      IF (D1.LT.EEWOR .AND. D2.LT.EEWOR)        GOTO 999
      X1     = XYHA(1,1)*COSI - XYHA(2,1)*SINI
      Y1     = XYHA(1,1)*SINI + XYHA(2,1)*COSI
      X2     = XYHA(1,2)*COSI - XYHA(2,2)*SINI
      Y2     = XYHA(1,2)*SINI + XYHA(2,2)*COSI
      H1     = XYHA(3,1)
      H2     = XYHA(3,2)
      IF (ABS(H1-H2) .LT. EEWOR)               GOTO 200
      DX1    = X1 - (X2*COSDA(K)-Y2*SINDA(K))
      DY1    = Y1 - (X2*SINDA(K)+Y2*COSDA(K))
      DH1    = H1 - H2
      DX2    = (X1*COSDA(K)-Y1*SINDA(K)) - X2
      DY2    = (X1*SINDA(K)+Y1*COSDA(K)) - Y2
      DH2    = H1 - H2
      AA     = DY1*DH2 - DY2*DH1
      BB     = DH1*DX2 - DH2*DX1
      CC     = DX1*DY2 - DX2*DY1
      S      = SQRT(AA*AA + BB*BB + CC*CC)
      IF (S .LT. EEWOR)                         GOTO 999
      AA     = AA / S
      BB     = BB / S
      CC     = CC / S
      DD     =-(AA*X1 + BB*Y1 + CC*H1)
      DO 130 I=1,NASTP(K)
        CG(J+KCGAF) = 0.
        CG(J+KCGAA) = AA
        CG(J+KCGBB) = BB
        CG(J+KCGCC) = CC
        CG(J+KCGDD) = DD
        AW     = AA*COSDA(K) - BB*SINDA(K)
        BB     = AA*SINDA(K) + BB*COSDA(K)
        AA     = AW
        CG(J+KCGNE) = 4.
        IF (D1.LT.EEWOR .OR. D2.LT.EEWOR)       CG(J+KCGNE) = 3.
        J      = J + LCGFAC
*
        CG(J+KCGAE) = ATREDG
        IF (I.EQ.1 .AND. IFULL.EQ.0)    CG(J+KCGAE) = 0.
        CG(J+KCGX1) = X2
        CG(J+KCGY1) = Y2
        CG(J+KCGZ1) = H2
        CG(J+KCGX2) = X1
        CG(J+KCGY2) = Y1
        CG(J+KCGZ2) = H1
        J      = J + LCGEDG
*
        IF (D1 .LT. EEWOR)              GOTO 110
        CG(J+KCGAE) = XYHA(4,1)
        CG(J+KCGX1) = X1
        CG(J+KCGY1) = Y1
        CG(J+KCGZ1) = H1
        IF (I.NE.NASTP(K) .OR.  IFULL.EQ.0)
     +                         CG(J+KCGX2) = X1*COSDA(K) - Y1*SINDA(K)
        IF (I.NE.NASTP(K) .OR.  IFULL.EQ.0)
     +                         CG(J+KCGY2) = X1*SINDA(K) + Y1*COSDA(K)
        IF (I.EQ.NASTP(K) .AND. IFULL.NE.0)
     +                         CG(J+KCGX2) = CG(LCGFAC+KCGX2)
        IF (I.EQ.NASTP(K) .AND. IFULL.NE.0)
     +                         CG(J+KCGY2) = CG(LCGFAC+KCGY2)
        CG(J+KCGZ2) = H1
        X1     = CG(J+KCGX2)
        Y1     = CG(J+KCGY2)
        J      = J + LCGEDG
*
  110   IF (D2 .LT. EEWOR)              GOTO 120
        CG(J+KCGAE) = XYHA(4,2)
        CG(J+KCGX2) = X2
        CG(J+KCGY2) = Y2
        CG(J+KCGZ2) = H2
        IF (I.NE.NASTP(K) .OR.  IFULL.EQ.0)
     +                         CG(J+KCGX1) = X2*COSDA(K) - Y2*SINDA(K)
        IF (I.NE.NASTP(K) .OR.  IFULL.EQ.0)
     +                         CG(J+KCGY1) = X2*SINDA(K) + Y2*COSDA(K)
        IF (I.EQ.NASTP(K) .AND. IFULL.NE.0)
     +                         CG(J+KCGX1) = CG(LCGFAC+KCGX1)
        IF (I.EQ.NASTP(K) .AND. IFULL.NE.0)
     +                         CG(J+KCGY1) = CG(LCGFAC+KCGY1)
        CG(J+KCGZ1) = H2
        X2     = CG(J+KCGX1)
        Y2     = CG(J+KCGY1)
        J      = J + LCGEDG
*
  120   CG(J+KCGAE) = ATREDG
        IF (I.EQ.NASTP(K) .AND. IFULL.EQ.0)     CG(J+KCGAE) = 0.
        CG(J+KCGX1) = X1
        CG(J+KCGY1) = Y1
        CG(J+KCGZ1) = H1
        CG(J+KCGX2) = X2
        CG(J+KCGY2) = Y2
        CG(J+KCGZ2) = H2
        J      = J + LCGEDG
  130   CONTINUE
      GOTO 999
*
**          M A K E   H O R I Z O N T A L   F A C E
*
  200 DD     = D1 - D2
      IF (ABS(DD) .LT. EEWOR)          GOTO 999
      K1     = K
      K2     = 3 - K
      CG(J+KCGAF) = 0.
      CG(J+KCGAA) = 0.
      CG(J+KCGBB) = 0.
      IF (DD .GT. 0)                    CG(J+KCGCC) = 1.
      IF (DD .LT. 0)                    CG(J+KCGCC) =-1.
      CG(J+KCGDD) =-XYHA(3,1) * CG(J+KCGCC)
      NEDGE  = 0.
      IF (D1 .GT. 0)                    NEDGE = NEDGE + NASTP(K1)
      IF (D2 .GT. 0)                    NEDGE = NEDGE + NASTP(K2)
      IF (IFULL .EQ. 0)                 NEDGE = NEDGE + 2
      CG(J+KCGNE) = NEDGE
      J      = J + LCGFAC
      IF (IFULL .NE. 0)                 GOTO 210
      CG(J+KCGAE) = 0.
      CG(J+KCGX1) = X2
      CG(J+KCGY1) = Y2
      CG(J+KCGZ1) = H2
      CG(J+KCGX2) = X1
      CG(J+KCGY2) = Y1
      CG(J+KCGZ2) = H1
      J      = J + LCGEDG
  210 IF (D1 .EQ. 0.)                    GOTO 230
      DO 220 I=1,NASTP(K1)
        CG(J+KCGAE) = XYHA(4,1)
        CG(J+KCGX1) = X1
        CG(J+KCGY1) = Y1
        CG(J+KCGZ1) = H1
        CG(J+KCGX2) = X1*COSDA(K1) - Y1*SINDA(K1)
        CG(J+KCGY2) = X1*SINDA(K1) + Y1*COSDA(K1)
        CG(J+KCGZ2) = H1
        X1     = CG(J+KCGX2)
        Y1     = CG(J+KCGY2)
        J      = J + LCGEDG
  220   CONTINUE
      IF (IFULL .EQ. 0)                  GOTO 230
      X1     = XYHA(1,1)*COSI - XYHA(2,1)*SINI
      Y1     = XYHA(1,1)*SINI + XYHA(2,1)*COSI
      CG(J-LCGEDG+KCGX2) = X1
      CG(J-LCGEDG+KCGY2) = Y1
*
  230 IF (D2 .EQ. 0.)                    GOTO 250
      DO 240 I=1,NASTP(K2)
        CG(J+KCGAE) = XYHA(4,2)
        CG(J+KCGX2) = X2
        CG(J+KCGY2) = Y2
        CG(J+KCGZ2) = H2
        CG(J+KCGX1) = X2*COSDA(K2) - Y2*SINDA(K2)
        CG(J+KCGY1) = X2*SINDA(K2) + Y2*COSDA(K2)
        CG(J+KCGZ1) = H2
        X2     = CG(J+KCGX1)
        Y2     = CG(J+KCGY1)
        J      = J + LCGEDG
  240   CONTINUE
      IF (IFULL .EQ. 0)                  GOTO 250
      X2     = XYHA(1,2)*COSI - XYHA(2,2)*SINI
      Y2     = XYHA(1,2)*SINI + XYHA(2,2)*COSI
      CG(J-LCGEDG+KCGX1) = X2
      CG(J-LCGEDG+KCGY1) = Y2
*
  250 IF (IFULL .NE. 0)                 GOTO 999
      CG(J+KCGAE) = 0.
      CG(J+KCGX1) = X1
      CG(J+KCGY1) = Y1
      CG(J+KCGZ1) = H1
      CG(J+KCGX2) = X2
      CG(J+KCGY2) = Y2
      CG(J+KCGZ2) = H2
      J      = J + LCGEDG
*
  999 LCG    = J
      RETURN
      END
