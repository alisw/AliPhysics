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
      SUBROUTINE CGHDFA(NFACE,XYZ,IFACE,IPFACE,DFACE)
************************************************************************
*                                                                      *
*     Name: CGHDFA                                                     *
*     Author: E. Chernyaev                       Date:    07.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Compute faces min-max                                  *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: NFACE - number of faces                                   *
*         XYZ(3,*) - node coordinates                                  *
*         IFACE(*) - face array                                        *
*        IPFACE(*) - pointers to faces beginings in IFACE              *
*                                                                      *
*     Output: DFACE(6,*) - faces min-max                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL      XYZ(3,*),DFACE(6,*)
*SG
      INTEGER      IFACE(*)
      INTEGER      IPFACE(*)
*
*  Error status flag for integer*2 variables :
*SG
      COMMON /CGOVER/ IHSERR
       IHSERR=0
*JS
*-
      DO 300 NF=1,NFACE
        J      = IPFACE(NF)
        NEDGE  = IFACE(J)
        J      = J + 1
        DFACE(1,NF) = 99999.
        DFACE(2,NF) = 99999.
        DFACE(3,NF) = 99999.
        DFACE(4,NF) =-99999.
        DFACE(5,NF) =-99999.
        DFACE(6,NF) =-99999.
        DO 200 NE=1,NEDGE
          DO 100 NN=1,2
            N      = IFACE(J)
***SG
*    Eliminating integer*2 variables, it's impossible to obtain
*    iface(j) < 0.
*JS
            IF(N.LT.0)THEN
               IHSERR=1
               RETURN
            ENDIF
*JS
***SG
            J      = J + 1
            IF(XYZ(1,N) .LT. DFACE(1,NF))       DFACE(1,NF) = XYZ(1,N)
            IF(XYZ(2,N) .LT. DFACE(2,NF))       DFACE(2,NF) = XYZ(2,N)
            IF(XYZ(3,N) .LT. DFACE(3,NF))       DFACE(3,NF) = XYZ(3,N)
            IF(XYZ(1,N) .GT. DFACE(4,NF))       DFACE(4,NF) = XYZ(1,N)
            IF(XYZ(2,N) .GT. DFACE(5,NF))       DFACE(5,NF) = XYZ(2,N)
            IF(XYZ(3,N) .GT. DFACE(6,NF))       DFACE(6,NF) = XYZ(3,N)
  100       CONTINUE
  200     CONTINUE
        DFACE(1,NF) =-DFACE(1,NF)
        DFACE(2,NF) =-DFACE(2,NF)
        DFACE(3,NF) = DFACE(4,NF)
        DFACE(4,NF) = DFACE(5,NF)
        DFACE(5,NF) = DFACE(6,NF)
        DFACE(6,NF) =-99999.
  300   CONTINUE
      RETURN
      END
