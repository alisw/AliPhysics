*
* $Id$
*
* $Log$
* Revision 1.1  1997/08/21 17:06:53  mclareni
* New routine from F.James
*
*
      SUBROUTINE RN3DIM(X,Y,Z,XLONG)
 
      DIMENSION RBF3(30),RBF2(20)
 
      SAVE IBF3,RBF3,IBF2,RBF2
 
      DATA IBF3 /30/, IBF2 /20/
 
    1 IF(IBF3 .GE. 30) THEN
       CALL RANLUX(RBF3,30)
       IBF3 = 1
      ENDIF
      A=RBF3(IBF3)-0.5
      B=RBF3(IBF3+1)-0.5
      C=RBF3(IBF3+2)-0.5
      IBF3=IBF3+3
      R=A**2+B**2+C**2
      IF(R .GT. 0.25) GO TO 1
      RINV=XLONG/SQRT(R)
      X=A*RINV
      Y=B*RINV
      Z=C*RINV
      RETURN
 
      ENTRY RN2DIM(X,Y,XLONG)
 
    2 IF(IBF2 .GE. 20) THEN
       CALL RANLUX(RBF2,20)
       IBF2 = 1
      ENDIF
      A=RBF2(IBF2)-0.5
      B=RBF2(IBF2+1)-0.5
      IBF2=IBF2+2
      R=A**2+B**2
      IF(R .GT. 0.25) GO TO 2
      RINV=XLONG/SQRT(R)
      X=A*RINV
      Y=B*RINV
      RETURN
      END
