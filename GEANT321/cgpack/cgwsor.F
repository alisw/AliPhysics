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
      SUBROUTINE CGWSOR(NWIRE,WIRE)
************************************************************************
*                                                                      *
*     Name: CGWSOR                                                     *
*     Author: E. Chernyaev                       Date:    16.04.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Prepare wire array                                     *
*                                                                      *
*     References: CGSHEL                                               *
*                                                                      *
*     Input:  NWIRE - number of edges                                  *
*     Output: WIRE(*,*)  - edges                                       *
*                                                                      *
*     Output:                                                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      WIRE(LCGEDG,*),SAVE(LCGEDG)
*-
      ERROR  = 1.5*EEWOR
      NE     = 0
      ATRINV = 0.
      CALL CGSHEL(ATRINV,ERROR,NWIRE,WIRE)
      IF (NWIRE .LE. 0)                 GOTO 999
*
**          D E L E T E   U N V I S I B L E    D O U B L E  E D G E S
**          A T R I B U T E   A N A L I S A T I O N
*
      J      = 1
  100 IF (J .GT. NWIRE)                 GOTO 999
      A1          = WIRE(KCGAE,J)
      SAVE(KCGX1) = WIRE(KCGX1,J)
      SAVE(KCGY1) = WIRE(KCGY1,J)
      SAVE(KCGZ1) = WIRE(KCGZ1,J)
      SAVE(KCGX2) = WIRE(KCGX2,J)
      SAVE(KCGY2) = WIRE(KCGY2,J)
      SAVE(KCGZ2) = WIRE(KCGZ2,J)
      J     = J + 1
  150 IF (J .GT. NWIRE)                                 GOTO 200
      IF (ABS(WIRE(KCGX1,J)-SAVE(KCGX1)) .GT. ERROR)    GOTO 200
      IF (ABS(WIRE(KCGY1,J)-SAVE(KCGY1)) .GT. ERROR)    GOTO 200
      IF (ABS(WIRE(KCGZ1,J)-SAVE(KCGZ1)) .GT. ERROR)    GOTO 200
      IF (ABS(WIRE(KCGX2,J)-SAVE(KCGX2)) .GT. ERROR)    GOTO 200
      IF (ABS(WIRE(KCGY2,J)-SAVE(KCGY2)) .GT. ERROR)    GOTO 200
      IF (ABS(WIRE(KCGZ2,J)-SAVE(KCGZ2)) .GT. ERROR)    GOTO 200
*           D O U B L E   W I R E
      A2     = WIRE(KCGAE,J)
      J      = J + 1
      IF (A1.LT.0. .AND. A2.LT.0.)                      GOTO 100
      IF (A1.LT.0.)                                     A1 =-A1 - 1.
      IF (A2.LT.0.)                                     A2 =-A2 - 1.
      IF (A2.GT.A1)                                     A1 = A2
      GOTO 150
*           I N S E R T   I N   W I R E (*)
  200 NE      = NE + 1
      IF (A1.LT.0.)                                     A1 =-A1 - 1.
      WIRE(KCGAE,NE) = A1
      WIRE(KCGX1,NE) = SAVE(KCGX1)
      WIRE(KCGY1,NE) = SAVE(KCGY1)
      WIRE(KCGZ1,NE) = SAVE(KCGZ1)
      WIRE(KCGX2,NE) = SAVE(KCGX2)
      WIRE(KCGY2,NE) = SAVE(KCGY2)
      WIRE(KCGZ2,NE) = SAVE(KCGZ2)
      GOTO 100
*
  999 NWIRE  = NE
      RETURN
      END
