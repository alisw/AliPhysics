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
      SUBROUTINE CGSUB(A,B,NMAX,C)
************************************************************************
*                                                                      *
*     Name: CGSUB                                                      *
*     Author: E. Chernyaev                       Date:    05.04.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Subrtarc one CG-objects from another                   *
*                                                                      *
*     References: CGBOOL                                               *
*                                                                      *
*     Input: A(*) - 1-st CG-object                                     *
*            B(*) - 2-nd CG-object                                     *
*            NMAX - max length C array                                 *
*                                                                      *
*     Output: C(*) - resulting CG-object                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL      A(*),B(*),C(*)
*-
      CALL CGBOOL(1,0,A,B,NMAX,C)
      RETURN
      END
