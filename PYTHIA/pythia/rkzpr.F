 
*==================================================================
 
      FUNCTION RKZPR(L,Q1,Q2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX RKZPR
      DIMENSION Q1(0:4),Q2(0:4)
      IF(IABS(L).NE.1) THEN
        WRITE(6,*) ' RKZPR: ERROR   L=',L
        STOP
      ENDIF
C...Introduce cutoff to check that R1 and R2 not zero.
      R1=DSQRT(MAX(1D-10,Q1(0)-Q1(1)))
      R2=DSQRT(MAX(1D-10,Q2(0)-Q2(1)))
      RKZPR=CMPLX(SNGL(Q1(2)),SNGL(Q1(3)))*R2/R1
     .     -CMPLX(SNGL(Q2(2)),SNGL(Q2(3)))*R1/R2
      IF(L.EQ.-1) RKZPR=-CONJG(RKZPR)
      END
