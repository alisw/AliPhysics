
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
      
      icrf=1
      irot=1
      
C---> Particle energies ---------
      P1S=P1X*P1X+P1Y*P1Y+P1Z*P1Z
      P2S=P2X*P2X+P2Y*P2Y+P2Z*P2Z
      E1=DSQRT(AM1*AM1+P1S)
      E2=DSQRT(AM2*AM2+P2S)
C---> Pair parameters -----------
      E12=E1+E2       ! Energy
      P12X=P1X+P2X    ! Px
      P12Y=P1Y+P2Y    ! Py
      P12Z=P1Z+P2Z    ! Pz
      P12S=P12X**2+P12Y**2+P12Z**2
      P12 =DSQRT(P12S)! Momentum
      V12 =P12/E12    ! Velocity
      CTH =P12Z/P12   ! cos(theta)
      STH =DSQRT(1.D0-CTH**2) !sin
      V12Z=V12*CTH    ! Longit. V
      GAMZ=1.D0/DSQRT(1.D0-V12Z**2)
C--      V12T=V12*STH    ! Transv. V in CMS (not needed)
      P12TS=P12X*P12X+P12Y*P12Y
      P12T=DSQRT(P12TS) !Pt
C===> Azimuthal rotation (Pt||X) ============
      IF(V12T.NE.0.D0) THEN
        CPHI=P12X/P12T  ! cos(phi)
        SPHI=P12Y/P12T  ! sin(phi)
        IF((irot.eq.1)) THEN 
          CALL ROT8(P1X,P1Y,SPHI,CPHI,P1X,P1Y)       
          CALL ROT8(P2X,P2Y,SPHI,CPHI,P2X,P2Y)
          CALL ROT8(X1,Y1,SPHI,CPHI,X1,Y1)       
          CALL ROT8(X2,Y2,SPHI,CPHI,X2,Y2) 
        END IF             
      ELSE ! Rotation impossible 
       CPHI=2.D0 ! to avoid
       SPHI=2.D0 ! using it !  
      END IF             
C===> Co-moving ref. frame       ============
      IF(icrf.gt.0) THEN
        CALL LTR8(P1Z,E1,V12Z,GAMZ,P1Z,E1a)
        CALL LTR8(P2Z,E2,V12Z,GAMZ,P2Z,E2a)
        P1S=P1X*P1X+P1Y*P1Y+P1Z*P1Z
        P2S=P2X*P2X+P2Y*P2Y+P2Z*P2Z
        E1=DSQRT(AM1*AM1+P1S)
        E2=DSQRT(AM2*AM2+P2S)
        CALL LTR8(Z1,T1,V12Z,GAMZ,Z1,T1)
        CALL LTR8(Z2,T2,V12Z,GAMZ,Z2,T2)
      END IF        
C===> Pair reference frame       ============
      P1=DSQRT(P1S)
      P2=DSQRT(P2S)
      E12=E1+E2
      P12X=P1X+P2X
      P12Y=P1Y+P2Y
      P12Z=P1Z+P2Z
      P12S=P12X**2+P12Y**2+P12Z**2
      P12 =DSQRT(P12S)
      AM12S=E12*E12-P12S
      AM12=DSQRT(AM12S)
      EPM=E12+AM12
      P112=P1X*P12X+P1Y*P12Y+P1Z*P12Z
      H1=(P112/EPM-E1)/AM12
      PPX=P1X+P12X*H1
      PPY=P1Y+P12Y*H1
      PPZ=P1Z+P12Z*H1
      EE=(E12*E1-P112)/AM12
      AKS=EE**2-AM1**2
      AK=DSQRT(AKS)
      CVK=(P12X*PPX+P12Y*PPY+P12Z*PPZ)/(P12*AK)
      V=P12/E12 
      V12T=P12T/SQRT(AM12S+P12TS) ! transverse velocity in LCMS   
C---> Coordinates -----------------------------
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


      RETURN
      END
C====       
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
      RPS=R
      RP=R*R
      RETURN
      END











 
