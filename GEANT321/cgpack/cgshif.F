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
      SUBROUTINE CGSHIF(DX,DY,DZ,CG)
************************************************************************
*                                                                      *
*     Name: CGSHIF                                                     *
*     Author: E. Chernyaev                       Date:    24.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Shift of  CG-object                                    *
*                                                                      *
*     References: CGAFFI                                               *
*                                                                      *
*     Input:     DX - X-shift                                          *
*                DY - Y-shift                                          *
*                DZ - Z-shift                                          *
*     Output: CG(*) - CG-object                                        *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL      CG(*),T(4,3)
*-
      T(1,1) = 1.
      T(2,1) = 0.
      T(3,1) = 0.
      T(4,1) = DX
*
      T(1,2) = 0.
      T(2,2) = 1.
      T(3,2) = 0.
      T(4,2) = DY
*
      T(1,3) = 0.
      T(2,3) = 0.
      T(3,3) = 1.
      T(4,3) = DZ
*
      CALL CGAFFI(T,CG)
      RETURN
      END
