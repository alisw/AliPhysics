 
      SUBROUTINE ltran12
C==>calculating particle momentum in PRF {EE,PPX,PPY,PPZ} from
C-  the momentum of the first particle {E1,P1X,P1Y,P1Z) in NRF
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, !4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_CONS/PI,PI2,SPI,DR,W
      COMMON/LEDWEIGHT/WEIF,WEI,WEIN,ITEST,IRANPOS
      
 
C   fm --> 1/GeV
c      write(*,*)'in LTRAN12 W', W
c      write(*,*)'in LTRAN12 p1x p1y p1z', P1X,P1Y,P1Z
c      write(*,*)'in LTRAN12 p2x p2y p2z', P2X,P2Y,P2Z
c       write(*,*)'in LTRAN12 x y z t', X,Y,Z,T
c      write(*,*)'in LTRAN12 x2 y2 z2', X2,Y2,Z2,T2

CMLV
       IF(IRANPOS.EQ.0)THEN      
      X1=X1*W
      Y1=Y1*W
      Z1=Z1*W
      T1=T1*W
      X2=X2*W
      Y2=Y2*W
      Z2=Z2*W
      T2=T2*W
C   calculating Ri, Pi and Ei
      R1=DSQRT(X1*X1+Y1*Y1+Z1*Z1)
      R2=DSQRT(X2*X2+Y2*Y2+Z2*Z2)
      ENDIF
CMLV

      P1S=P1X*P1X+P1Y*P1Y+P1Z*P1Z
      P2S=P2X*P2X+P2Y*P2Y+P2Z*P2Z
      P1=DSQRT(P1S)
      P2=DSQRT(P2S)
      E1=DSQRT(AM1*AM1+P1S)
      E2=DSQRT(AM2*AM2+P2S)
C-----------------------------------------------------------------------
      E12=E1+E2
      P12X=P1X+P2X
      P12Y=P1Y+P2Y
      P12Z=P1Z+P2Z
      P12S=P12X**2+P12Y**2+P12Z**2
      AM12=DSQRT(E12**2-P12S)
      EPM=E12+AM12
      P12=DSQRT(P12S)
      P112=P1X*P12X+P1Y*P12Y+P1Z*P12Z
      H1=(P112/EPM-E1)/AM12
      PPX=P1X+P12X*H1
      PPY=P1Y+P12Y*H1
      PPZ=P1Z+P12Z*H1
      EE=(E12*E1-P112)/AM12
      AKS=EE**2-AM1**2
      AK=DSQRT(AKS)
CW      WRITE(6,38)'AK ',AK,'K ',PPX,PPY,PPZ,EE
38    FORMAT(A7,E11.4,A7,4E11.4)
      RETURN
      END
 