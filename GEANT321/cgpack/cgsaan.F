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
      SUBROUTINE CGSAAN(A1,A2,NA,NB,IREP)
************************************************************************
*                                                                      *
*     Name: CGSAAN                                                     *
*     Author: E. Chernyaev                       Date:    01.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Set parameters of rotation                             *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: A1 - intitial angle                                       *
*            A2 - end angle                                            *
*            NA - 1-st number of steps                                 *
*            NB - 2-nd number of steps                                 *
*                                                                      *
*     Output: IREP - reply (-1 if error in parameters)                 *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgcaan.inc"
      SAVE IFIRST
      DATA      IFIRST/1/
*-
      IF (IFIRST .NE. 0)                RAD = ATAN(1.0D+00)*4./180.
      IFIRST = 0
      IREP   =-1
      IF (NA .LE. 0)                    GOTO 999
      IF (NB .LE. 0)                    GOTO 999
      IF (A1 .GE. A2)                   GOTO 999
      DA     = A2 - A1
      IF (DA .GT. 360.)                 GOTO 999
      IF (DA .EQ. 360.)                 IFULL = 1
      IF (DA .LT. 360.)                 IFULL = 0
      NASTP(1) = NA
      NASTP(2) = NB
      COSI     = COS(A1*RAD)
      SINI     = SIN(A1*RAD)
      COSDA(1) = COS(DA*RAD/NA)
      SINDA(1) = SIN(DA*RAD/NA)
      COSDA(2) = COS(DA*RAD/NB)
      SINDA(2) = SIN(DA*RAD/NB)
      IREP   = 0
  999 RETURN
      END
