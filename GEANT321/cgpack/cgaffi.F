*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:41  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGAFFI(TT,CG)
************************************************************************
*                                                                      *
*     Name: CGAFFI                                                     *
*     Author: E. Chernyaev                       Date:    24.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Affine transformation of CG-object                     *
*                                                                      *
*     References: CGTSTR                                               *
*                                                                      *
*     Input:  T(4,3) - transformation matrix                           *
*     Output:  CG(*) - CG-object                                       *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      TT(4,3),CG(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION T(4,3),X1,Y1,Z1,X2,Y2,Z2,
     +                 A,B,C,Q,S,XGRAV,YGRAV,ZGRAV
#endif
#if defined(CERNLIB_SINGLE)
      REAL      T(4,3)
#endif
*-
      TDEL   = EEWOR / 10.
      CALL CGTSTR(CG,IREP)
      IF (IREP .LT. 0)          GOTO 999
      DO 20 J=1,3
        DO 10 I=1,4
          T(I,J) = TT(I,J)
          IF (ABS(TT(I,J))    .LT. TDEL)        T(I,J) = 0.
          IF (ABS(1.-TT(I,J)) .LT. TDEL)        T(I,J) = 1.
          IF (ABS(1.+TT(I,J)) .LT. TDEL)        T(I,J) =-1.
   10     CONTINUE
   20   CONTINUE
*
**          T R A N S F E R   C O O R D I N A T E S
*
      NFACE  = CG(KCGNF)
      IF (NFACE .EQ. 0)         GOTO 999
      JCG    = LCGHEA
      DO 200 NF=1,NFACE
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC
        DO 100 NE=1,NEDGE
          X1     = CG(JCG+KCGX1)
          Y1     = CG(JCG+KCGY1)
          Z1     = CG(JCG+KCGZ1)
          X2     = CG(JCG+KCGX2)
          Y2     = CG(JCG+KCGY2)
          Z2     = CG(JCG+KCGZ2)
          CG(JCG+KCGX1) = T(1,1)*X1 + T(2,1)*Y1 + T(3,1)*Z1 + T(4,1)
          CG(JCG+KCGY1) = T(1,2)*X1 + T(2,2)*Y1 + T(3,2)*Z1 + T(4,2)
          CG(JCG+KCGZ1) = T(1,3)*X1 + T(2,3)*Y1 + T(3,3)*Z1 + T(4,3)
          CG(JCG+KCGX2) = T(1,1)*X2 + T(2,1)*Y2 + T(3,1)*Z2 + T(4,1)
          CG(JCG+KCGY2) = T(1,2)*X2 + T(2,2)*Y2 + T(3,2)*Z2 + T(4,2)
          CG(JCG+KCGZ2) = T(1,3)*X2 + T(2,3)*Y2 + T(3,3)*Z2 + T(4,3)
          JCG    = JCG + LCGEDG
  100     CONTINUE
  200   CONTINUE
*
**          S E T   N E W   N O R M A L E S
*
      JCG    = LCGHEA
      DO 400 NF=1,NFACE
        JCGFAC = JCG
        NEDGE  = CG(JCG+KCGNE)
        A      = (T(2,2)*T(3,3) - T(3,2)*T(2,3))*CG(JCG+KCGAA) +
     +           (T(3,2)*T(1,3) - T(1,2)*T(3,3))*CG(JCG+KCGBB) +
     +           (T(1,2)*T(2,3) - T(2,2)*T(1,3))*CG(JCG+KCGCC)
        B      = (T(2,3)*T(3,1) - T(3,3)*T(2,1))*CG(JCG+KCGAA) +
     +           (T(3,3)*T(1,1) - T(1,3)*T(3,1))*CG(JCG+KCGBB) +
     +           (T(1,3)*T(2,1) - T(2,3)*T(1,1))*CG(JCG+KCGCC)
        C      = (T(2,1)*T(3,2) - T(3,1)*T(2,2))*CG(JCG+KCGAA) +
     +           (T(3,1)*T(1,2) - T(1,1)*T(3,2))*CG(JCG+KCGBB) +
     +           (T(1,1)*T(2,2) - T(2,1)*T(1,2))*CG(JCG+KCGCC)
*        IF (ABS(A) .LT. EEWOR/10.)   A=0.
*        IF (ABS(B) .LT. EEWOR/10.)   B=0.
*        IF (ABS(C) .LT. EEWOR/10.)   C=0.
        Q      = A*A + B*B + C*C
        S      = 1.
        IF (Q.GT.1+TDEL .OR. Q.LT.1-TDEL)       S = SQRT(Q)
        IF(S.LT.TDEL)S=TDEL
        XGRAV  = 0.
        YGRAV  = 0.
        ZGRAV  = 0.
        JCG    = JCG + LCGFAC
        DO 300 NE=1,NEDGE
          XGRAV  = XGRAV + CG(JCG+KCGX1) + CG(JCG+KCGX2)
          YGRAV  = YGRAV + CG(JCG+KCGY1) + CG(JCG+KCGY2)
          ZGRAV  = ZGRAV + CG(JCG+KCGZ1) + CG(JCG+KCGZ2)
          JCG    = JCG + LCGEDG
  300     CONTINUE
        XGRAV  = XGRAV / (2*NEDGE)
        YGRAV  = YGRAV / (2*NEDGE)
        ZGRAV  = ZGRAV / (2*NEDGE)
        CG(JCGFAC+KCGAA) = A / S
        CG(JCGFAC+KCGBB) = B / S
        CG(JCGFAC+KCGCC) = C / S
        CG(JCGFAC+KCGDD) =-(A*XGRAV + B*YGRAV + C*ZGRAV) / S
  400   CONTINUE
*
  999 RETURN
      END
