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
      SUBROUTINE CGHINS(CG,H,ISHAPE)
************************************************************************
*                                                                      *
*     Name: CGHINS                                                     *
*     Author: E. Chernyaev                       Date:    04.08.88     *
*     Revised by: S.Giani                        Date:    24.04.91     *
*                                                                      *
*     Function: Add next CG-object to HIDE-structure                   *
*                                                                      *
*     References: CGHSIZ, CGHMOV                                       *
*                                                                      *
*     Input: CG(*) - CG-object                                         *
*             H(*) - HIDE-structure                                    *
*                                                                      *
*     Output: none                                                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cghpar.inc"
*SG
#include "geant321/gcgobj.inc"
*SG
      REAL CG(*),H(*)
*-          P A R A M E T E R S   C H E C K
      IF (H(KHSIZE)  .LE. 0.)   GOTO 999
      CALL CGTSTR(CG,IREP)
*      IF (IREP .LT. 0) PRINT *, 'CGHINS: CG-structure test failed'
      IF (IREP .LT. 0)          GOTO 999
*
      NHSIZE = H(KHSIZE)
      NT     = H(KHNT)
      NXYZ   = H(KHNXYZ)
      NFACE  = H(KHNFAC)
      NEDGE  = H(KHNEDG)
      CALL CGHSIZ(CG,NT,NVF,NVE,ISHAPE)
      IF (NVF .EQ. 0)           GOTO 999
      NX     = 3*(NXYZ + NVE*2)
      NF     = (NFACE + NEDGE*2 + NVF + NVE*2)/I4SIZE + 1
      NP     = (NVE*2*2)/I4SIZE + 1
      NTOTAL = LHHEAD + NX + NF + NP
*           T E S T   S P A C E
      IF (NTOTAL .LT. NHSIZE)   GOTO 100
*SG
*   Setting new number of words needed
*
      NCLAS1=NTOTAL+1
      KCGST=-8
*SG
      H(KHNXYZ) = NXYZ + NVE*2
      H(KHNFAC) = NFACE + NVF
      H(KHNEDG) = NEDGE + NVE
      GOTO 999
*
  100 JXYZ1   = LHHEAD + NXYZ*3 + 1
      JPOINT = JXYZ1 + NVE*2*3
      JFACE  = 1
*SG
*   Resetting flag for size evaluation of Wire Structure
*
      KCGST=0
*SG
      CALL CGHMOV(CG,H,NVF,NVE,H(JXYZ1),H(JPOINT),H(JFACE))
*
  999 RETURN
      END
