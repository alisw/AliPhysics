*CMZ :          17/07/98  15.44.33  by  Federico Carminati
*-- Author :
C*********************************************************************

      SUBROUTINE LUROBO(THE,PHI,BEX,BEY,BEZ)

C...Purpose: to perform rotations and boosts.
      IMPLICIT DOUBLE PRECISION(D)
*KEEP,LUJETS.
      COMMON /LUJETS/ N,K(200000,5),P(200000,5),V(200000,5)
      SAVE /LUJETS/
*KEEP,LUDAT1.
      COMMON /LUDAT1/ MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /LUDAT1/
*KEND.
      DIMENSION ROT(3,3),PR(3),VR(3),DP(4),DV(4)

C...Find range of rotation/boost. Convert boost to double precision.
      IMIN=1
      IF(MSTU(1).GT.0) IMIN=MSTU(1)
      IMAX=N
      IF(MSTU(2).GT.0) IMAX=MSTU(2)
      DBX=BEX
      DBY=BEY
      DBZ=BEZ
      GOTO 110

C...Entry for specific range and double precision boost.
      ENTRY LUDBRB(IMI,IMA,THE,PHI,DBEX,DBEY,DBEZ)
      IMIN=IMI
      IF(IMIN.LE.0) IMIN=1
      IMAX=IMA
      IF(IMAX.LE.0) IMAX=N
      DBX=DBEX
      DBY=DBEY
      DBZ=DBEZ

C...Optional resetting of V (when not set before.)
      IF(MSTU(33).NE.0) THEN
        DO 100 I=MIN(IMIN,MSTU(4)),MIN(IMAX,MSTU(4))
        DO 100 J=1,5
  100   V(I,J)=0.
        MSTU(33)=0
      ENDIF

C...Check range of rotation/boost.
  110 IF(IMIN.GT.MSTU(4).OR.IMAX.GT.MSTU(4)) THEN
        CALL LUERRM(11,'(LUROBO:) range outside LUJETS memory')
        RETURN
      ENDIF

C...Rotate, typically from z axis to direction (theta,phi).
      IF(THE**2+PHI**2.GT.1E-20) THEN
        ROT(1,1)=COS(THE)*COS(PHI)
        ROT(1,2)=-SIN(PHI)
        ROT(1,3)=SIN(THE)*COS(PHI)
        ROT(2,1)=COS(THE)*SIN(PHI)
        ROT(2,2)=COS(PHI)
        ROT(2,3)=SIN(THE)*SIN(PHI)
        ROT(3,1)=-SIN(THE)
        ROT(3,2)=0.
        ROT(3,3)=COS(THE)
        DO 140 I=IMIN,IMAX
        IF(K(I,1).LE.0) GOTO 140
        DO 120 J=1,3
        PR(J)=P(I,J)
  120   VR(J)=V(I,J)
        DO 130 J=1,3
        P(I,J)=ROT(J,1)*PR(1)+ROT(J,2)*PR(2)+ROT(J,3)*PR(3)
  130   V(I,J)=ROT(J,1)*VR(1)+ROT(J,2)*VR(2)+ROT(J,3)*VR(3)
  140   CONTINUE
      ENDIF

C...Boost, typically from rest to momentum/energy=beta.
      IF(DBX**2+DBY**2+DBZ**2.GT.1E-20) THEN
        DB=SQRT(DBX**2+DBY**2+DBZ**2)
        IF(DB.GT.0.99999999D0) THEN
C...Rescale boost vector if too close to unity.
          CALL LUERRM(3,'(LUROBO:) boost vector too large')
          DBX=DBX*(0.99999999D0/DB)
          DBY=DBY*(0.99999999D0/DB)
          DBZ=DBZ*(0.99999999D0/DB)
          DB=0.99999999D0
        ENDIF
        DGA=1D0/SQRT(1D0-DB**2)
        DO 160 I=IMIN,IMAX
        IF(K(I,1).LE.0) GOTO 160
        DO 150 J=1,4
        DP(J)=P(I,J)
  150   DV(J)=V(I,J)
        DBP=DBX*DP(1)+DBY*DP(2)+DBZ*DP(3)
        DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP(4))
        P(I,1)=DP(1)+DGABP*DBX
        P(I,2)=DP(2)+DGABP*DBY
        P(I,3)=DP(3)+DGABP*DBZ
        P(I,4)=DGA*(DP(4)+DBP)
        DBV=DBX*DV(1)+DBY*DV(2)+DBZ*DV(3)
        DGABV=DGA*(DGA*DBV/(1D0+DGA)+DV(4))
        V(I,1)=DV(1)+DGABV*DBX
        V(I,2)=DV(2)+DGABV*DBY
        V(I,3)=DV(3)+DGABV*DBZ
        V(I,4)=DGA*(DV(4)+DBV)
  160   CONTINUE
      ENDIF

      RETURN
      END
