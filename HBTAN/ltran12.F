       SUBROUTINE ltran12
C

c      SUBROUTINE TRANS(icrf,irot)
C==>  TRANSformation to the Co-moving frame  (icrf>0) and
C     Rotation to the system where (Pt || X),(irot=1).
C-FSI ***************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  ! momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. 
     1                X2,Y2,Z2,T2,R2  ! points in NRF
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! momenta in PRF
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_POC/AMN,AM1,AM2,CN,C1,C2,AC1,AC2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM
      COMMON/FSI_CVK/V,CVK
C-FSI ***************************************************
      COMMON /PAIR/P12T,V12Z,GAMZ,V12T,CPHI,SPHI
 
 
C   calculating Ri, Pi and Ei
      R1=DSQRT(X1*X1+Y1*Y1+Z1*Z1)
      R2=DSQRT(X2*X2+Y2*Y2+Z2*Z2)
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

      RETURN
      END
C====       
C==== ===============================================================
C==== ===============================================================
C==== ===============================================================

      subroutine BoostToPrf()
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_CVK/V,CVK
      COMMON/FSI_MOM/P1X,P1Y,P1Z,E1,P1,  !part. momenta in NRF
     1               P2X,P2Y,P2Z,E2,P2
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS,
     1               X,Y,Z,T,RP,RPS
      COMMON/FSI_COOR/X1,Y1,Z1,T1,R1, ! 4-coord. of emis. points in NRF
     1                X2,Y2,Z2,T2,R2
      COMMON/FSI_P12/P12X,P12Y,P12Z,E12,P12,AM12,EPM

      XS=X1-X2
      YS=Y1-Y2
      ZS=Z1-Z2
      TS=T1-T2
      RS12=XS*P12X+YS*P12Y+ZS*P12Z
      H1=(RS12/EPM-TS)/AM12
      X=XS+P12X*H1
      Y=YS+P12Y*H1
      Z=ZS+P12Z*H1
      T=(E12*TS-RS12)/AM12
      RPS=X*X+Y*Y+Z*Z
      RP=DSQRT(RPS)
CW      WRITE(6,38)'RP ',RP,'X ',X,Y,Z,T
38    FORMAT(A7,E11.4,A7,4E11.4)
 
      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12
      return 
      end
C==== ===============================================================
C==== ===============================================================

      SUBROUTINE LTR8(Z,T,BETA,GAMMA,ZT,TT)
C===> Lorentz Transf. of Z(Pz) and T(E) to moving ref. frame.(REAL*8) 
CInp: Z,T-Zcoord,Time before tr., BETA,GAMMA- velocity, Lor.fact.
COut: ZT,TT- "    "   after transformation.
C==== ===============================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      ZH=GAMMA*(Z-BETA*T)     
      TT=GAMMA*(T-BETA*Z)
      ZT=ZH
      RETURN
      END
C====       
      
      SUBROUTINE LTR4(Z,T,BETA,GAMMA,ZT,TT)
C===> Lorentz Transf. of Z(Pz) and T(E) to moving ref. frame.(real*4) 
CInp: Z,T-Zcoord,Time before tr., BETA,GAMMA- velocity, Lor.fact.
COut: ZT,TT- "    "   after transformation.
C==== ===============================================================
      ZH=GAMMA*(Z-BETA*T)     
      TT=GAMMA*(T-BETA*Z)
      ZT=ZH
      RETURN
      END
C====           
      SUBROUTINE ROT8(X,Y,SF,CF,XR,YR)
C===> Rotation with the angle f. (REAL*8)
CInp: X,Y-coord. before rotation; SF=sin(f), CF=cos(f),
COut: XR,YR - coordinates after rotation.
C==== =================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      XH=X*CF+Y*SF !Y              
      YR=Y*CF-X*SF !   _-X'   
      XR=XH        ! _- f      
      RETURN       !------>
      END          !     X

      SUBROUTINE SETPDIST(R)
C=====Just sets distance between particles
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FSI_PRF/PPX,PPY,PPZ,AK,AKS, ! momenta in PRF
     1               X,Y,Z,T,RP,RPS
      RP=R
      RPS=R*R
      RETURN
      END
