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
      SUBROUTINE CGSPHE(R,KA,KB,LCG,CG)
************************************************************************
*                                                                      *
*     Name: CGSPHE                                                     *
*     Author: E. Chernyaev                       Date:    24.01.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Create CG-object for SPHERE (BALL)                     *
*                                                                      *
*     References: CGELLI                                               *
*                                                                      *
*     Input:  R  - sphere radius                                       *
*             KA - number of latitude step                             *
*             KB - number of longitude step                            *
*            LCG - max-size of CG-object                               *
*                                                                      *
*     Output: CG - CG-object                                           *
*                  CG(1) - length of CG-object                         *
*                        = 0 if error in parameters                    *
*                        < 0 if no space                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL      CG(*)
*-
      CALL CGELLI(R,R,R,KA,KB,LCG,CG)
      RETURN
      END
