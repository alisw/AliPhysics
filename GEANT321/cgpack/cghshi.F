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
      SUBROUTINE CGHSHI(NSIZE,NFACE,NEDGE,IHIDE,IFACE,IPFACE)
************************************************************************
*                                                                      *
*     Name: CGHSHI                                                     *
*     Author: E. Chernyaev                       Date:    07.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Shift face array and set pointers to faces             *
*                                                                      *
*     References:                                                      *
*                                                                      *
*     Input : NSIZE - size of IHIDE array (in INTEGER*2 words)         *
*     Output: NFACE - number of faces before and after shift           *
*             NEDGE - number of edges before and after shift           *
*          IHIDE(*) - old face array                                   *
*          IFACE(*) - new face array                                   *
*         IPFACE(*) - pointers to faces beginings in IFACE             *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
*SG
      INTEGER  IHIDE(*),IFACE(*)
      INTEGER  IPFACE(*)
*SG
*-
      NF     = 0
      NE     = 0
      JHIDE  = NSIZE - (NFACE+NEDGE*2) + 1
      JFACE  = 1
      DO 200 NNF=1,NFACE
        IE     = 0
        IEDGE  = IHIDE(JHIDE)
        JHIDE  = JHIDE + 1
        JWF    = JFACE + 1
        DO 100 NNE=1,IEDGE
          N1     = IHIDE(JHIDE)
          N2     = IHIDE(JHIDE+1)
          JHIDE  = JHIDE + 2
          IF (N1 .EQ. N2)       GOTO 100
          IFACE(JWF)   = N1
          IFACE(JWF+1) = N2
          JWF    = JWF + 2
          IE     = IE + 1
  100     CONTINUE
        IF (IE .EQ. 0)          GOTO 200
        NF     = NF + 1
        NE     = NE + IE
        IFACE(JFACE) = IE
        IPFACE(NF)   = JFACE
        JFACE  = JWF
  200   CONTINUE
      NFACE  = NF
      NEDGE  = NE
      RETURN
      END
