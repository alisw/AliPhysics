      FUNCTION VDOT (X,Y,N)
C
C CERN PROGLIB# F121    VDOT            .VERSION KERNFOR  1.0   710701
C ORIG. 01/07/71
C
      DIMENSION X(*),Y(*)
C
      XX= 0.
      IF (N.LE.0)  GO TO 100
      DO 9 I= 1,N
    9 XX= XX + X(I)*Y(I)
C
  100 VDOT= XX
      RETURN
      END
