*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:38  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.25  by  S.Giani
*-- Author :
      FUNCTION GPSCIN(EBEG,EEND,ECUR,COEFF)
C.
C.    ******************************************************************
C.    *                                                                *
C.    ******************************************************************
C.
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
      REAL COEFF(4)
      PARAMETER (ONE=1)
C.
C.    ------------------------------------------------------------------
C.
      YBEG = EBEG/ECUR
      YEND = EEND/ECUR
      YBEGIN = ONE/YBEG
      YENDIN = ONE/YEND
*
      EINV = ONE/ECUR
      C1 = COEFF(1)*EINV
      C2 = COEFF(2)*EINV**2
      C3 = COEFF(3)*EINV**3
      C4 = COEFF(4)*EINV**4
*
      CONS1 = 0.5*(C1+C3)
      CONS2 = C2+C4
      CONS3 = 0.5*CONS2
      CONS4 = 0.5*C3
      CONS5 = C4/3
*
*     VALBEG = CONS1*F1(YBEG)+CONS3*F2(YBEG)+
*    +         YBEGIN*(CONS2+YBEGIN*(CONS4+YBEGIN*CONS5))
*     VALEND = CONS1*F1(YEND)+CONS3*F2(YEND)+
*    +         YENDIN*(CONS2+YENDIN*(CONS4+YENDIN*CONS5))
*
*     GPSCIN = VALEND-VALBEG
      GPSCIN =
     +      CONS1*LOG(ABS((YEND**2-1)*YBEG**2/((YBEG**2-1)*YEND**2)))+
     +      CONS2*(YENDIN-YBEGIN)+
     +      CONS3*LOG(ABS((1-YEND)*(1+YBEG)/((1-YBEG)*(1+YEND))))+
     +      CONS4*(YENDIN**2-YBEGIN**2)+CONS5*(YENDIN**3-YBEGIN**3)
*
      END
