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
      SUBROUTINE CGSCOP(IFIRST,CG,RMIN,RMAX)
************************************************************************
*                                                                      *
*     Name: CGSCOP                                                     *
*     Author: E. Chernyaev                       Date:    30.01.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Compute CUBE scope of set of CG-objects                *
*                                                                      *
*     References: CGMNMX                                               *
*                                                                      *
*     Input:  IFIRST - flag (1 - if begin of the set)                  *
*             CG(*) - CG-object                                        *
*                                                                      *
*     Output: RMIN(3) - min coordinates of the CUBE                     *
*             RMAX(3) - max coordinates of the CUBE                     *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL      CG(*),RMIN(3),RMAX(3)
      REAL      RMID(3),RCMIN(3),RCMAX(3)
*-
      CALL CGMNMX(CG,RMIN,RMAX)
      IF (IFIRST .EQ. 0)                GOTO 100
      RCMIN(1) = RMIN(1)
      RCMIN(2) = RMIN(2)
      RCMIN(3) = RMIN(3)
      RCMAX(1) = RMAX(1)
      RCMAX(2) = RMAX(2)
      RCMAX(3) = RMAX(3)
  100 IF (RMIN(1) .GT. RMAX(1))         GOTO 999
      RRR = RCMAX(1) - RCMIN(1)
      DO 200 I=1,3
        IF (RCMIN(I) .GT. RMIN(I))      RCMIN(I) = RMIN(I)
        IF (RCMAX(I) .LT. RMAX(I))      RCMAX(I) = RMAX(I)
        RMID(I) = (RCMIN(I) + RCMAX(I))/2.
        IF (RRR .LT. RCMAX(I)-RCMIN(I)) RRR = RCMAX(I) - RCMIN(I)
  200   CONTINUE
      DO 300 I=1,3
        RCMIN(I) = RMID(I) - RRR/2.
        RCMAX(I) = RMID(I) + RRR/2.
  300   CONTINUE
*
  999 RMIN(1)  = RCMIN(1)
      RMIN(2)  = RCMIN(2)
      RMIN(3)  = RCMIN(3)
      RMAX(1)  = RCMAX(1)
      RMAX(2)  = RCMAX(2)
      RMAX(3)  = RCMAX(3)
      RETURN
      END
