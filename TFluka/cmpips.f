*$ CREATE CMPIPS.FOR
*COPY CMPIPS
*                                                                      *
*=== Cmpips ===========================================================*
*                                                                      *
      SUBROUTINE CMPIPS ( NALLD0, STEPID, STEPTT, STEP )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 2005-2006:                Alfredo Ferrari          *
*     All rights reserved                                              *
*                                                                      *
*                                                                      *
*     CoMpute Primary Ionization PoSitions:                            *
*                                                                      *
*     Created  on  11 november 2005  by       Alfredo Ferrari          *
*                                               INFN - Milan           *
*                                                                      *
*     Last change  on  05-oct-06     by  Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*     Input variables:                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(ALLDLT)'
      INCLUDE '(TRACKR)'
*
      LOGICAL LBGSTP
      DIMENSION CUMTTR (0:MXTRCK), RNDVEC (MXALLD), INDVEC (MXALLD)
*
      SAVE NSTART, CUMTTR
*
      LBGSTP = STEPTT .LT. AZRZRZ
      if (ntrack .eq. 0) THEN
*         WRITE(6,*) "Warning ntrack = 0", NALLD0, STEPID, STEPTT, STEP
         RETURN
      ENDIF
*  +-------------------------------------------------------------------*
*  |  Beginning of a step:
      IF ( LBGSTP ) THEN
         NSTART = 1
         CUMTTR (0) = ZERZER
         DO 1000 I = 1, NTRACK
            CUMTTR (I) = CUMTTR (I-1) + TTRACK (I)
 1000    CONTINUE
         CRVCRR = STEP / CUMTTR (NTRACK)
      END IF
*  |
*  +-------------------------------------------------------------------*
      NRNGEN = MIN ( NALLDL, MXALLD ) - NALLD0
      IF ( NRNGEN .LE. 0 ) RETURN
      CALL FLRNLP ( RNDVEC, NRNGEN )
*  The previous line can be substituted by:
*     DO 1400 I = 1, NRNGEN
*        RNDVEC (I) = FLRNDM (RNDPOI)
*1400 CONTINUE
      CALL RORDIN ( RNDVEC, INDVEC, NRNGEN )
*  +-------------------------------------------------------------------*
*  |  Loop on primary electrons:
      DO 5000 I = 1, NRNGEN
         TTRACR = ( STEPTT + RNDVEC (I) * STEPID ) / CRVCRR
         NSTAR0 = NSTART
         DO 3000 J = NSTAR0, NTRACK
            IF ( TTRACR .LT. CUMTTR (J) ) THEN
               NSTART = J
               REDUC  = ( CUMTTR (J) - TTRACR ) / TTRACK (J)
               GO TO 4000
            END IF
 3000    CONTINUE
         CALL FLABRT ( 'CMPIPS', '3000-CONTINUE' )
 4000    CONTINUE
         IF ( REDUC .LT. ZERZER .OR. REDUC .GT. ONEONE )
     &      CALL FLABRT ( 'CMPIPS', 'INVALID REDUC' )
         K = I + NALLD0
         XALLDL (K) = REDUC * ( XTRACK (NSTART) - XTRACK (NSTART-1) )
     &              + XTRACK (NSTART-1)
         YALLDL (K) = REDUC * ( YTRACK (NSTART) - YTRACK (NSTART-1) )
     &              + YTRACK (NSTART-1)
         ZALLDL (K) = REDUC * ( ZTRACK (NSTART) - ZTRACK (NSTART-1) )
     &              + ZTRACK (NSTART-1)
 5000 CONTINUE
*  |
*  +-------------------------------------------------------------------*
      RETURN
*=== End of subroutine Cmpips =========================================*
      END

