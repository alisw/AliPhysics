*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGHCRE(NT,IFWIND,W1,W2,N,H)
************************************************************************
*                                                                      *
*     Name: CGHCRE                                                     *
*     Author: E. Chernyaev                       Date:    05.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Initiate HIDE-structure                                *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input :   NT - number of transformation to screen system         *
*           IFWIND - window flag                                       *
*            W1(2) - left lower window corner                          *
*            W2(2) - rigth upper window corner                         *
*                N - max size of H array                               *
*             H(*) - HIDE-structure                                    *
*                                                                      *
*     Output: none                                                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cghpar.inc"
#include "geant321/cgctra.inc"
      REAL H(*),W1(2),W2(2)
*-
      H(KHSIZE) = 0.
      IF (N .LT. LHHEAD)                GOTO 999
      DO 100 I=1,LHHEAD
        H(I)    = 0.
  100   CONTINUE
      IF (NT.LE.0 .OR. NT.GT.NTMAX)     GOTO 999
      H(KHSIZE) = N
      H(KHNT)   = NT
      H(KHIFWI) = IFWIND
      IF (IFWIND .LE. 0)                GOTO 999
      H(KHIFWI+1) = W1(1)
      H(KHIFWI+2) = W1(2)
      H(KHIFWI+3) = W2(1)
      H(KHIFWI+4) = W2(2)
  999 RETURN
      END
