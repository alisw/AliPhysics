*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      FUNCTION ULANGL(X,Y)

C...Purpose: to reconstruct an angle from given x and y coordinates.
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.

      ULANGL=0.
      R=SQRT(X**2+Y**2)
      IF(R.LT.1E-20) RETURN
      IF(ABS(X)/R.LT.0.8) THEN
        ULANGL=SIGN(ACOS(X/R),Y)
      ELSE
        ULANGL=ASIN(Y/R)
        IF(X.LT.0..AND.ULANGL.GE.0.) THEN
          ULANGL=PARU(1)-ULANGL
        ELSEIF(X.LT.0.) THEN
          ULANGL=-PARU(1)-ULANGL
        ENDIF
      ENDIF

      RETURN
      END
