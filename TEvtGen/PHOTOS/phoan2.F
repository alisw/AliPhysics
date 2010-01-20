      FUNCTION PHOAN2(X,Y)
C.----------------------------------------------------------------------
C.
C.    PHOTOS:   PHOton radiation in decays calculation of ANgle '2'
C.
C.    Purpose:  Calculate angle from X and Y
C.
C.    Input Parameters:  X, Y
C.
C.    Output Parameter:  Function value
C.
C.    Author(s):  S. Jadach                       Created at:  01/01/89
C.                B. van Eijk                     Last Update: 02/01/90
C.
C.----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION PHOAN2
      REAL*8 X,Y
      REAL*8 PI,TWOPI
      COMMON/PHPICO/PI,TWOPI
      IF (ABS(Y).LT.ABS(X)) THEN
        PHOAN2=ATAN(ABS(Y/X))
        IF (X.LE.0.D0) PHOAN2=PI-PHOAN2
      ELSE
        PHOAN2=ACOS(X/SQRT(X**2+Y**2))
      ENDIF
      RETURN
      END
