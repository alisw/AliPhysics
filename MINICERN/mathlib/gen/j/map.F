*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:51  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
       SUBROUTINE MAP(F,XL,XU,XD)
C      WRITE(6,)A TABLE AND PLOT OF  THE FUNCTION F
      DIMENSION Y(512),LINE(103)
      EXTERNAL F
      WRITE(6,500)
 500  FORMAT('Q')
#if defined(CERNLIB_IBM)||defined(CERNLIB_CDC)
      I=LOCF(F)-LOCF(LINE)
#endif
#if defined(CERNLIB_CDC)
      WRITE(6,1)LINE(I-2)
    1 FORMAT('1',40('.'),'   TABLE AND PLOT OF THE FUNCTION  ',A7)
#endif
#if defined(CERNLIB_IBM)
      I=LINE(I+1)/4-LOCF(LINE)
      WRITE(6,1)LINE(I+2),LINE(I+3)
    1 FORMAT('1',40('.'),'   TABLE AND PLOT OF THE FUNCTION  ',A4,A3)
#endif
#if defined(CERNLIB_VAX)
      WRITE(6,1)
    1 FORMAT('1',40('.'),'   TABLE AND PLOT OF THE USER FUNCTION')
#endif
      IF  (XU .LE. XL)  GO TO 20
      IF  (XD .LE. 0.)  GO TO 30
      IF  ( (XU-XL)/XD .LT. 2.)  GO TO 40
       DO 2 I=1,512
       X=XL+(I-1)*XD
       IF(X.GT.XU)GO TO 3
    2 Y(I)=F(X)
       I=513
    3 NX=I-1
       YXU=F(XU)
       YL=Y(1)
       YU=YL
       DO 4 I=2,NX
       Y1=Y(I)
       YL=MIN(Y1,YL)
    4 YU=MAX(Y1,YU)
       NY=100
       YD=(YU-YL)/NY
       WRITE(6,5)XL,YL,XU,YU,XD,YD,NX,NY
    5 FORMAT('0',14X,'ABSCISSA PARAMETERS=',48X,'ORDINATE PARAMETERS='/
     +  15X,'LOWER LIMIT...',E12.3,42X,'LOWER LIMIT...',E12.3/
     +  15X,'UPPER LIMIT...',E12.3,42X,'UPPER LIMIT...',E12.3/
     +  15X,'STEP SIZE.....',E12.3,42X,'BIN SIZE......',E12.3/
     +  15X,'POINT NUMBER..',  I12,42X,'BIN NUMBER....',I12/)
      WRITE(6,6)(I,I=10,100,10)
    6 FORMAT('0POINT  ABSCISSA  ORDINATE    0',9I10,I11)
      YP=YL*YU
      IF(YP.GE.0)GO TO 101
      K=INT(-YL/YD)+2
  101 CONTINUE
      WRITE(6,9)
    9 FORMAT(29X,'O',10('+',9('-')),'+O')
       YS=0.0
      CALL VBLANK(LINE,103)
       DO 11 I=1,NX
       N=I-1
       X=XL+N*XD
       Y1=Y(I)
       YS=YS+Y1
       J=INT((Y1-YL)/YD)+2
      IF(YP.LT.0.) CALL VFILL(LINE(K),1,1H.)
      CALL VFILL(LINE(1),1,1HI)
      LINE(103)=LINE(1)
       IF(MOD(N,10).NE.0)GO TO 12
      CALL VFILL(LINE(1),1,1H+)
      LINE(103)=LINE(1)
   12 CALL VFILL(LINE(J),1,1H*)
       WRITE(6,13)N,X,Y1,LINE
   13 FORMAT(I6,1X,2E10.3,2X,103A1)
   11 CALL VBLANK(LINE(J),1)
      WRITE(6,9)
       V=2.0*YS-Y(1)-Y(NX)
       V=0.5*XD*V
       V=V+0.5*(XU-X)*(Y(NX)+YXU)
       WRITE(6,16)V
   16 FORMAT('0','ESTIMATE OF INTEGRAL...',E10.3)
      WRITE(6,17)
   17 FORMAT('0',130('='))
       RETURN
 20   WRITE(6,200)
      STOP 1
 30   WRITE(6,300)
      STOP 1
 40   WRITE(6,400)
      STOP 1
 200  FORMAT('0UPPER LIMIT <= LOWER LIMIT - MAP STOPS')
 300  FORMAT('0DISTANCE BETWEEN POINTS <= 0 - MAP STOPS')
 400  FORMAT('0DISTANCE BETWEEN POINTS > INTERVAL/2 - MAP STOPS')
       END
