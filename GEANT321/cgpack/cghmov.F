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
      SUBROUTINE CGHMOV(CG,H,NVF,NVE,XYZ,IP,IFACE)
************************************************************************
*                                                                      *
*     Name: CGHMOV                                                     *
*     Author: E. Chernyaev                       Date:    04.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Move visible faces from CG-object to HIDE-structure    *
*                                                                      *
*     References: CGHREN                                               *
*                                                                      *
*     Input :   CG - CG-object                                         *
*                H - HIDE-structure                                    *
*              NVF - number of visible faces                           *
*              NVE - number of visible edges                           *
*         XYZ(3,*) - array for new nodes coordinates                   *
*          IP(2,*) - work array for node renumbering                   *
*         IFACE(*) - faces array                                       *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cghpar.inc"
      REAL      CG(*),H(*),XYZ(6,*)
*SG
      INTEGER   IP(2,*),IFACE(*)
*SG
*-
      JFACE  = H(KHSIZE)*I4SIZE-((H(KHNFAC)+NVF)+(H(KHNEDG)+NVE)*2)+1
      JF     = 0
      JE     = 0
      NN     = 0
      JCG    = LCGHEA
      NFACE  = CG(KCGNF)
      DO 400 NF=1,NFACE
        NEDGE = CG(JCG+KCGNE)
        IF (NEDGE .LT. 0)       GOTO 300
        IFACE(JFACE+JF) = NEDGE
        DO 100 I=1,NEDGE*2
          IFACE(JFACE+JF+I) = NN + I
  100     CONTINUE
        DO 200 I=1,NEDGE
          J     = JCG + LCGFAC + (I-1)*LCGEDG
          XYZ(1,JE+I) = CG(J+KCGX1)
          XYZ(2,JE+I) = CG(J+KCGY1)
          XYZ(3,JE+I) = CG(J+KCGZ1)
          XYZ(4,JE+I) = CG(J+KCGX2)
          XYZ(5,JE+I) = CG(J+KCGY2)
          XYZ(6,JE+I) = CG(J+KCGZ2)
  200     CONTINUE
        NN    = NN + NEDGE*2
        JF    = JF + NEDGE*2 + 1
        JE    = JE + NEDGE
  300   IF (NEDGE .LT. 0)       NEDGE =-NEDGE
        CG(JCG+KCGNE) = NEDGE
        JCG   = JCG + LCGFAC + NEDGE*LCGEDG
  400   CONTINUE
      NT     = H(KHNT)
      NOLD   = H(KHNXYZ)
      NNEW   = NVE*2
      CALL CGHREN(NT,NOLD,NNEW,NVF,XYZ,IP,IFACE(JFACE),NTOTAL)
      H(KHNXYZ) = NTOTAL
      H(KHNFAC) = H(KHNFAC) + NVF
      H(KHNEDG) = H(KHNEDG) + NVE
*
  999 RETURN
      END
