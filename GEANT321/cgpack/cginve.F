*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGINVE(CG)
************************************************************************
*                                                                      *
*     Name: CGINVE                                                     *
*     Author: E. Chernyaev                       Date:    20.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Invert edge direction                                  *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: CG(*) - CG-object                                         *
*     Output:                                                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      CG(*)
*-
      NFACE  = CG(KCGNF)
      JCG    = LCGHEA
      DO 200 NF=1,NFACE
        CG(JCG+KCGAA) =-CG(JCG+KCGAA)
        CG(JCG+KCGBB) =-CG(JCG+KCGBB)
        CG(JCG+KCGCC) =-CG(JCG+KCGCC)
        CG(JCG+KCGDD) =-CG(JCG+KCGDD)
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC
        DO 100 NE=1,NEDGE
          X      = CG(JCG+KCGX1)
          Y      = CG(JCG+KCGY1)
          Z      = CG(JCG+KCGZ1)
          CG(JCG+KCGX1) = CG(JCG+KCGX2)
          CG(JCG+KCGY1) = CG(JCG+KCGY2)
          CG(JCG+KCGZ1) = CG(JCG+KCGZ2)
          CG(JCG+KCGX2) = X
          CG(JCG+KCGY2) = Y
          CG(JCG+KCGZ2) = Z
          JCG    = JCG + LCGEDG
  100     CONTINUE
  200   CONTINUE
      RETURN
      END
