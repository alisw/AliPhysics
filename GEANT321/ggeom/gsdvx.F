*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GSDVX (NAME,MOTHER,NDIV,IAXIS,STEP,C0,NUMED
     +                  ,NDVMX)
C.
C.    *******************************************************
C.    *                                                     *
C.    * DIVIDES MOTHER INTO 'NDIV' DIVISIONS OF SIZE 'STEP' *
C.    * ALONG AXIS 'IAXIS' STARTING AT COORDINATE 'C0'.     *
C.    *  NDIV OR STEP CAN BE SET NEGATIVE OR ZERO, IN WHICH *
C.    * CASE THEY ARE RECOMPUTED FROM THE MOTHER'S SIZE.    *
C.    *                                                     *
C.    *    ==>Called by : 'USER'                            *
C.    *    Authors   F.Bruyant, B.Stoer  ********           *
C.    *                                                     *
C.    *******************************************************
C.
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
      CHARACTER*4 NAME,MOTHER
C.    -------------------------------------------------------
*
      IF (NDIV.GT.0) THEN
         IF (STEP.GT.0.) THEN
            CALL GSDVT2 (NAME, MOTHER, STEP, IAXIS, C0, NUMED, -NDIV)
         ELSE
            CALL GSDVN2 (NAME, MOTHER, NDIV, IAXIS, C0, NUMED)
         ENDIF
      ELSE
         IF (STEP.GT.0.) THEN
            CALL GSDVT2 (NAME, MOTHER, STEP, IAXIS, C0, NUMED, NDVMX)
         ELSE
            WRITE (CHMAIL, 1001)
            CALL GMAIL (0, 0)
            IEORUN = 1
         ENDIF
      ENDIF
*
 1001 FORMAT (' GSDVX : Call with NDIV and STEP .le. 0')
      END
