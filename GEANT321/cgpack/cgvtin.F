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
      SUBROUTINE CGVTIN(TMIN,TMAX,IVIS)
************************************************************************
*                                                                      *
*     Name: CGVTIN                                                     *
*     Author: E. Chernyaev                       Date:    18.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Insert new portion in hidden portions list             *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: TMIN - portion begin                                      *
*            TMAX - portion end                                        *
*                                                                      *
*     Output: IVIS = -1 if edge fully unvisible                        *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgcedg.inc"
*-
      IF (TMAX-TMIN .LT. TDEL)          GOTO 999
      KJ     = 0
      IF (NPART .EQ. 0)                 GOTO 210
*           C H E C K   C O V E R E D   P O R T I O N   OF LINE AGAINST
*           P R E V I O U S   C O V E R E D   P O R T I O N S
      DO 100 I = 1,NPART
        IF (TMIN.LE.TSTRT(I) .AND. KJ.EQ.0)     KJ = I
        IF (TMIN.GT.TEND(I)+TDEL .OR. TMAX.LT.TSTRT(I)-TDEL)    GOTO 100
*           T H E R E   I S   C O V E R E D   PORTIONS  INTERSECTIONS
        IF (TMIN.GE.TSTRT(I)-TDEL .AND. TMAX.LE.TEND(I)+TDEL)   GOTO 999
        IF (TMIN .LT. TSTRT(I))         TSTRT(I) = TMIN
        IF (TMAX .GT. TEND(I))          TEND(I)  = TMAX
        GOTO 300
  100   CONTINUE
*           N O   C O V E R E D   P O R T I O N S   INTERSECTION.
*           C R E A T E   N E W   P O R T I O N
      IF (KJ .EQ .0)                    GOTO 210
      DO 200 I = NPART,KJ,-1
        TSTRT(I+1) = TSTRT(I)
        TEND(I+1)  = TEND(I)
  200   CONTINUE
  210 NPART  = NPART + 1
      IF (KJ .EQ. 0)                    KJ = NPART
      TSTRT(KJ) = TMIN
      TEND(KJ)  = TMAX
*           C O V E R E D   P O R T I O N S   RECONSTRUCTION (IF NEED)
  300 IF (NPART .LE. 1)                 GOTO 600
      J      = 1
      DO 500 I = 2,NPART
        IF (TSTRT(I)-TEND(J) .LT. TDEL) GOTO 490
        J = J+1
        TSTRT(J) = TSTRT(I)
        TEND(J)  = TEND (I)
        GOTO 500
  490   TEND(J)  = MAX (TEND(J),TEND(I))
  500   CONTINUE
      NPART  = J
*           C H E C K   I F   L I N E   I S   F U L L Y   C O V E R E D
  600 T0     = 0.
      T1     = 1.
      IF (TSTRT(1) .LE. TDEL)           T0 = TEND(1)
      IF (TEND(NPART) .GE. 1.-TDEL)     T1 = TSTRT(NPART)
      IF (T0 .GE. T1)                   IVIS = -1
*
  999 RETURN
      END
