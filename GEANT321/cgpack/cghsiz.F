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
      SUBROUTINE CGHSIZ(CG,NT,NVF,NVE,ISHAPE)
************************************************************************
*                                                                      *
*     Name: CGHSIZ                                                     *
*     Author: E. Chernyaev                       Date:    02.08.88     *
*     Revised by: S.Giani                        Date:    26.04.91     *
*                                                                      *
*     Function: Compute number of visible FACES and EDGE               *
*                                                                      *
*     References: CGFVIS                                               *
*                                                                      *
*     Input:  CG - CG-object                                           *
*             NT - number for transformation to screen coordinates     *
*                                                                      *
*     Output: NVF - number of visible faces                            *
*             NVE - number of visible edges                            *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/gcspee.inc"
      REAL      CG(*)
*-
**SG
      NTFLAG = 1
      NTIM   = 0
**SG
      NVF    = 0
      NVE    = 0
      JCG    = LCGHEA
      NFACE  = CG(KCGNF)
      DO 100 NF=1,NFACE
        NEDGE = CG(JCG+KCGNE)
        CALL CGFARE(NT,CG(JCG+1),IVIS,ISHAPE)
        IF (IVIS .LT. 0)        CG(JCG+KCGNE) =-CG(JCG+KCGNE)
        JCG   = JCG + LCGFAC + NEDGE*LCGEDG
        IF (IVIS .LT. 0)        GOTO 100
        NVF   = NVF + 1
        NVE   = NVE + NEDGE
  100   CONTINUE
*SG
      NTFLAG = 0
*SG
  999 RETURN
      END
