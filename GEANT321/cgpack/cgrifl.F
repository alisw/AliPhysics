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
*
      SUBROUTINE CGRIFL(TT,CG)
*
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      TT(4,3),CG(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION T(4,3),DET1,DET2,DET3,DETTT
#endif
#if defined(CERNLIB_SINGLE)
      REAL      T(4,3)
#endif
*-
*      TDEL   = EEWOR / 10.
      CALL CGTSTR(CG,IREP)
      IF (IREP .LT. 0)          GOTO 999
      DO 20 J=1,3
        DO 10 I=1,4
          T(I,J) = TT(I,J)
*          IF (ABS(TT(I,J))    .LT. TDEL)        T(I,J) = 0.
*          IF (ABS(1.-TT(I,J)) .LT. TDEL)        T(I,J) = 1.
*          IF (ABS(1.+TT(I,J)) .LT. TDEL)        T(I,J) =-1.
   10     CONTINUE
   20   CONTINUE
*
      DET1=T(1,1)*(T(2,2)*T(3,3)-T(2,3)*T(3,2))
      DET2=T(1,2)*(T(2,1)*T(3,3)-T(2,3)*T(3,1))
      DET3=T(1,3)*(T(2,1)*T(3,2)-T(3,1)*T(2,2))
      DETTT=DET1-DET2+DET3
      IF(DETTT.LT.0.)THEN
*           S E T   I N V E R S E   O R D E R   O F   E D G E S
      JCG=LCGHEA
      NFACE  = CG(KCGNF)
      DO 500 NF=1,NFACE
        CG(JCG + KCGAA) =-CG(JCG + KCGAA)
        CG(JCG + KCGBB) =-CG(JCG + KCGBB)
        CG(JCG + KCGCC) =-CG(JCG + KCGCC)
        CG(JCG + KCGDD) =-CG(JCG + KCGDD)
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC
        DO 400 NE=1,NEDGE
          X      = CG(JCG+KCGX1)
          Y      = CG(JCG+KCGY1)
          Z      = CG(JCG+KCGZ1)
          CG(JCG + KCGX1) = CG(JCG + KCGX2)
          CG(JCG + KCGY1) = CG(JCG + KCGY2)
          CG(JCG + KCGZ1) = CG(JCG + KCGZ2)
          CG(JCG + KCGX2) = X
          CG(JCG + KCGY2) = Y
          CG(JCG + KCGZ2) = Z
          JCG    = JCG + LCGEDG
  400     CONTINUE
  500   CONTINUE
        ENDIF
  999   RETURN
        END
 
 
