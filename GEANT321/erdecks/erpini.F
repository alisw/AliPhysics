*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 15:37:35  mclareni
* Add geane321 source directories
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.49  by  S.Giani
*-- Author :
      SUBROUTINE ERPINI
*
C.    ******************************************************************
C.    *                                                                *
C.    *       Initializes the error propagation                        *
C.    *                                                                *
C.    *    ==>Called by : ERTRAK                                       *
C.    *       Author    E.Nagy  *********                              *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gckine.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcunit.inc"
#include "geant321/erwork.inc"
#include "geant321/ertrio.inc"
*
      DIMENSION DUMM(15)
*
* *** Starting point
*
      DO 11 I = 1,3
         XI(I) = VECT(I)
         PPI(I) = VECT(7)*VECT(3+I)
   11 CONTINUE
*
* *** Initialize tracking
*
      MVAR   = 0
      IFLAG  =-1
      ITRAN  = 0
      CALL TRPROP (XI, PPI, HI, DUM, DUM, HI, CHTR, DUM, DUMM,
     &            MVAR, IFLAG, ITRAN, IERR)
      IF (IERR.NE.0) THEN
         WRITE (LOUT, 777) IERR
         RETURN
      ENDIF
      DO 13 I = 1,15
         EF(I) = EI(I)
   13 CONTINUE
*
  777 FORMAT(//4X,' **** Error in subr. TRPROP',I4,'  called by'
     F     ,' subr. ERPINI'//)
*
       END
