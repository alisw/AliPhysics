* One step extrapolation routines in Fortran:
* EXTRAP_ONESTEP_HELIX, EXTRAP_ONESTEP_HELIX3, EXTRAP_ONESTEP_RUNGEKUTTA.
* Taken respectively from Geant (gtrak) routines:
* GHELIX, GHELIX3, GRKUTA.
* Everything in double precision,
* in order that the track fit with Minuit is converging.
* Modifications from Geant are indicated with "Cmodif".

Cmodif: SUBROUTINE GHELIX (CHARGE, STEP, VECT, VOUT) changed into:
      SUBROUTINE EXTRAP_ONESTEP_HELIX (CHARGE, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Performs the tracking of one step in a magnetic field         *
C.    *  The trajectory is assumed to be a helix in a constant field   *
C.    *  taken at the mid point of the step.                           *
C.    *  Parameters:                                                   *
C.    *   input                                                        *
C.    *     STEP =arc length of the step asked                         *
C.    *     VECT =input vector (position,direction cos and momentum)   *
C.    *     CHARGE=  electric charge of the particle                   *
C.    *   output                                                       *
C.    *     VOUT = same as VECT after completion of the step           *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Author    M.Hansroul  *********                          *
C.    *       Modified  S.Egli, S.V.Levonian                           *
C.    *       Modified  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.

Cmodif: everything in double precision
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION      VECT(7),VOUT(7)
      DIMENSION      XYZ(3),H(4),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      IF (CHARGE.EQ.0.)         GO TO 10
      XYZ(1)    = VECT(IX) + 0.5 * STEP * VECT(IPX)
      XYZ(2)    = VECT(IY) + 0.5 * STEP * VECT(IPY)
      XYZ(3)    = VECT(IZ) + 0.5 * STEP * VECT(IPZ)
C
Cmodif: CALL GUFLD (XYZ, H) changed into:
      CALL GUFLD_DOUBLE (XYZ, H)
 
      H2XY = H(1)**2 + H(2)**2
      H(4) = H(3)**2 + H2XY
      IF (H(4).LE.1.E-12)       GO TO 10
      IF (H2XY.LE.1.E-12*H(4))  THEN
Cmodif: CALL GHELX3 (CHARGE*H(3), STEP, VECT, VOUT) changed into:
         CALL EXTRAP_ONESTEP_HELIX3 (CHARGE*H(3), STEP, VECT, VOUT)
         GO TO 999
      ENDIF
      H(4) = SQRT(H(4))
      H(1) = H(1)/H(4)
      H(2) = H(2)/H(4)
      H(3) = H(3)/H(4)
      H(4) = H(4)*EC
*
      HXP(1) = H(2)*VECT(IPZ) - H(3)*VECT(IPY)
      HXP(2) = H(3)*VECT(IPX) - H(1)*VECT(IPZ)
      HXP(3) = H(1)*VECT(IPY) - H(2)*VECT(IPX)
 
      HP = H(1)*VECT(IPX) + H(2)*VECT(IPY) + H(3)*VECT(IPZ)
*
      RHO = -CHARGE*H(4)/VECT(IPP)
      TET = RHO * STEP
      IF (ABS(TET).GT.0.15)     THEN
         SINT = SIN(TET)
         SINTT = (SINT/TET)
         TSINT = (TET-SINT)/TET
         COS1T = 2.*(SIN(0.5*TET))**2/TET
      ELSE
         TSINT = SIXTH*TET**2
         SINTT = (1. - TSINT)
         SINT = TET*SINTT
         COS1T = 0.5*TET
      ENDIF
*
      F1 = STEP * SINTT
      F2 = STEP * COS1T
      F3 = STEP * TSINT * HP
      F4 = -TET*COS1T
      F5 = SINT
      F6 = TET * COS1T * HP
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1) + F3*H(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2) + F3*H(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F2*HXP(3) + F3*H(3))
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1) + F6*H(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2) + F6*H(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F5*HXP(3) + F6*H(3))
 
      GO TO 999
 
   10 CONTINUE
      DO 20 I   = 1,3
         VOUT(I) = VECT(I) + STEP * VECT(I+3)
         VOUT(I+3) = VECT(I+3)
   20 CONTINUE
C
  999 END

Cmodif: SUBROUTINE GHELX3 (FIELD, STEP, VECT, VOUT) changed into:
      SUBROUTINE EXTRAP_ONESTEP_HELIX3 (FIELD, STEP, VECT, VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Tracking routine in a constant field oriented            *
C.    *       along axis 3                                             *
C.    *       Tracking is performed with a conventional                *
C.    *       helix step method                                        *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Authors    R.Brun, M.Hansroul  *********                 *
C     *       Rewritten  V.Perevoztchikov
C.    *                                                                *
C.    ******************************************************************
C.

Cmodif: everything in double precision
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION      VECT(7),VOUT(7),HXP(3)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6,IPP=7)
      PARAMETER      (SIXTH = 1./6.)
      PARAMETER      (EC=2.9979251E-4)
C.
C.    ------------------------------------------------------------------
C.
C       units are kgauss,centimeters,gev/c
C
      VOUT(IPP) = VECT(IPP)
      H4 = FIELD * EC
*
      HXP(1) = - VECT(IPY)
      HXP(2) = + VECT(IPX)
 
      HP = VECT(IPZ)
*
      RHO = -H4/VECT(IPP)
      TET = RHO * STEP
      IF (ABS(TET).GT.0.15)     THEN
         SINT = SIN(TET)
         SINTT = (SINT/TET)
         TSINT = (TET-SINT)/TET
         COS1T = 2.*(SIN(0.5*TET))**2/TET
      ELSE
         TSINT = SIXTH*TET**2
         SINTT = (1. - TSINT)
         SINT = TET*SINTT
         COS1T = 0.5*TET
      ENDIF
*
      F1 = STEP * SINTT
      F2 = STEP * COS1T
      F3 = STEP * TSINT * HP
      F4 = -TET*COS1T
      F5 = SINT
      F6 = TET * COS1T * HP
 
      VOUT(IX) = VECT(IX) + (F1*VECT(IPX) + F2*HXP(1))
      VOUT(IY) = VECT(IY) + (F1*VECT(IPY) + F2*HXP(2))
      VOUT(IZ) = VECT(IZ) + (F1*VECT(IPZ) + F3)
 
      VOUT(IPX) = VECT(IPX) + (F4*VECT(IPX) + F5*HXP(1))
      VOUT(IPY) = VECT(IPY) + (F4*VECT(IPY) + F5*HXP(2))
      VOUT(IPZ) = VECT(IPZ) + (F4*VECT(IPZ) + F6)
 
C
  999 END

Cmodif: SUBROUTINE GRKUTA (CHARGE,STEP,VECT,VOUT) changed into:
      SUBROUTINE EXTRAP_ONESTEP_RUNGEKUTTA (CHARGE,STEP,VECT,VOUT)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Runge-Kutta method for tracking a particle through a magnetic *
C.    *  field. Uses Nystroem algorithm (See Handbook Nat. Bur. of     *
C.    *  Standards, procedure 25.5.20)                                 *
C.    *                                                                *
C.    *  Input parameters                                              *
C.    *       CHARGE    Particle charge                                *
C.    *       STEP      Step size                                      *
C.    *       VECT      Initial co-ords,direction cosines,momentum     *
C.    *  Output parameters                                             *
C.    *       VOUT      Output co-ords,direction cosines,momentum      *
C.    *  User routine called                                           *
C.    *       CALL GUFLD(X,F)                                          *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUSWIM                               *
C.    *       Authors    R.Brun, M.Hansroul  *********                 *
C.    *                  V.Perevoztchikov (CUT STEP implementation)    *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
Cmodif: no condition from CERNLIB_SINGLE for double precision
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
Cmodif: REAL changed into DOUBLE PRECISION in following 2 lines
      DOUBLE PRECISION CHARGE, STEP, VECT(*), VOUT(*), F(4)
      DOUBLE PRECISION XYZT(3), XYZ(3), X, Y, Z, XT, YT, ZT
      DIMENSION SECXS(4),SECYS(4),SECZS(4),HXP(3)
      EQUIVALENCE (X,XYZ(1)),(Y,XYZ(2)),(Z,XYZ(3)),
     +            (XT,XYZT(1)),(YT,XYZT(2)),(ZT,XYZT(3))
*
      PARAMETER (MAXIT = 1992, MAXCUT = 11)
      PARAMETER (EC=2.9979251D-4,DLT=1D-4,DLT32=DLT/32)
      PARAMETER (ZERO=0, ONE=1, TWO=2, THREE=3)
      PARAMETER (THIRD=ONE/THREE, HALF=ONE/TWO)
      PARAMETER (PISQUA=.986960440109D+01)
      PARAMETER      (IX=1,IY=2,IZ=3,IPX=4,IPY=5,IPZ=6)
*.
*.    ------------------------------------------------------------------
*.
*             This constant is for units CM,GEV/C and KGAUSS
*
      ITER = 0
      NCUT = 0
      DO 10 J=1,7
         VOUT(J)=VECT(J)
   10 CONTINUE
      PINV   = EC * CHARGE / VECT(7)
      TL = 0.
      H      = STEP
*
*
   20 REST  = STEP-TL
      IF (ABS(H).GT.ABS(REST)) H = REST
Cmodif: CALL GUFLD(VOUT,F) changed into:
      CALL GUFLD_DOUBLE(VOUT,F)
*
*             Start of integration
*
      X      = VOUT(1)
      Y      = VOUT(2)
      Z      = VOUT(3)
      A      = VOUT(4)
      B      = VOUT(5)
      C      = VOUT(6)
*
      H2     = HALF * H
      H4     = HALF * H2
      PH     = PINV * H
      PH2    = HALF * PH
      SECXS(1) = (B * F(3) - C * F(2)) * PH2
      SECYS(1) = (C * F(1) - A * F(3)) * PH2
      SECZS(1) = (A * F(2) - B * F(1)) * PH2
      ANG2 = (SECXS(1)**2 + SECYS(1)**2 + SECZS(1)**2)
      IF (ANG2.GT.PISQUA) GO TO 40
      DXT    = H2 * A + H4 * SECXS(1)
      DYT    = H2 * B + H4 * SECYS(1)
      DZT    = H2 * C + H4 * SECZS(1)
      XT     = X + DXT
      YT     = Y + DYT
      ZT     = Z + DZT
*
*              Second intermediate point
*
      EST = ABS(DXT)+ABS(DYT)+ABS(DZT)
      IF (EST.GT.H) GO TO 30
 
Cmodif: CALL GUFLD(XYZT,F) changed into:
      CALL GUFLD_DOUBLE(XYZT,F)
      AT     = A + SECXS(1)
      BT     = B + SECYS(1)
      CT     = C + SECZS(1)
*
      SECXS(2) = (BT * F(3) - CT * F(2)) * PH2
      SECYS(2) = (CT * F(1) - AT * F(3)) * PH2
      SECZS(2) = (AT * F(2) - BT * F(1)) * PH2
      AT     = A + SECXS(2)
      BT     = B + SECYS(2)
      CT     = C + SECZS(2)
      SECXS(3) = (BT * F(3) - CT * F(2)) * PH2
      SECYS(3) = (CT * F(1) - AT * F(3)) * PH2
      SECZS(3) = (AT * F(2) - BT * F(1)) * PH2
      DXT    = H * (A + SECXS(3))
      DYT    = H * (B + SECYS(3))
      DZT    = H * (C + SECZS(3))
      XT     = X + DXT
      YT     = Y + DYT
      ZT     = Z + DZT
      AT     = A + TWO*SECXS(3)
      BT     = B + TWO*SECYS(3)
      CT     = C + TWO*SECZS(3)
*
      EST = ABS(DXT)+ABS(DYT)+ABS(DZT)
      IF (EST.GT.2.*ABS(H)) GO TO 30
 
Cmodif: CALL GUFLD(XYZT,F) changed into:
      CALL GUFLD_DOUBLE(XYZT,F)
*
      Z      = Z + (C + (SECZS(1) + SECZS(2) + SECZS(3)) * THIRD) * H
      Y      = Y + (B + (SECYS(1) + SECYS(2) + SECYS(3)) * THIRD) * H
      X      = X + (A + (SECXS(1) + SECXS(2) + SECXS(3)) * THIRD) * H
*
      SECXS(4) = (BT*F(3) - CT*F(2))* PH2
      SECYS(4) = (CT*F(1) - AT*F(3))* PH2
      SECZS(4) = (AT*F(2) - BT*F(1))* PH2
      A      = A+(SECXS(1)+SECXS(4)+TWO * (SECXS(2)+SECXS(3))) * THIRD
      B      = B+(SECYS(1)+SECYS(4)+TWO * (SECYS(2)+SECYS(3))) * THIRD
      C      = C+(SECZS(1)+SECZS(4)+TWO * (SECZS(2)+SECZS(3))) * THIRD
*
      EST    = ABS(SECXS(1)+SECXS(4) - (SECXS(2)+SECXS(3)))
     ++        ABS(SECYS(1)+SECYS(4) - (SECYS(2)+SECYS(3)))
     ++        ABS(SECZS(1)+SECZS(4) - (SECZS(2)+SECZS(3)))
*
      IF (EST.GT.DLT .AND. ABS(H).GT.1.E-4) GO TO 30
      ITER = ITER + 1
      NCUT = 0
*               If too many iterations, go to HELIX
      IF (ITER.GT.MAXIT) GO TO 40
*
      TL = TL + H
      IF (EST.LT.(DLT32)) THEN
         H = H*TWO
      ENDIF
      CBA    = ONE/ SQRT(A*A + B*B + C*C)
      VOUT(1) = X
      VOUT(2) = Y
      VOUT(3) = Z
      VOUT(4) = CBA*A
      VOUT(5) = CBA*B
      VOUT(6) = CBA*C
      REST = STEP - TL
      IF (STEP.LT.0.) REST = -REST
      IF (REST .GT. 1.E-5*ABS(STEP)) GO TO 20
*
      GO TO 999
*
**              CUT STEP
   30 NCUT = NCUT + 1
*               If too many cuts , go to HELIX
      IF (NCUT.GT.MAXCUT)       GO TO 40
      H = H*HALF
      GO TO 20
*
**              ANGLE TOO BIG, USE HELIX
   40 F1  = F(1)
      F2  = F(2)
      F3  = F(3)
      F4  = SQRT(F1**2+F2**2+F3**2)
      RHO = -F4*PINV
      TET = RHO * STEP
      IF(TET.NE.0.) THEN
         HNORM = ONE/F4
         F1 = F1*HNORM
         F2 = F2*HNORM
         F3 = F3*HNORM
*
         HXP(1) = F2*VECT(IPZ) - F3*VECT(IPY)
         HXP(2) = F3*VECT(IPX) - F1*VECT(IPZ)
         HXP(3) = F1*VECT(IPY) - F2*VECT(IPX)
 
         HP = F1*VECT(IPX) + F2*VECT(IPY) + F3*VECT(IPZ)
*
         RHO1 = ONE/RHO
         SINT = SIN(TET)
         COST = TWO*SIN(HALF*TET)**2
*
         G1 = SINT*RHO1
         G2 = COST*RHO1
         G3 = (TET-SINT) * HP*RHO1
         G4 = -COST
         G5 = SINT
         G6 = COST * HP
 
         VOUT(IX) = VECT(IX) + (G1*VECT(IPX) + G2*HXP(1) + G3*F1)
         VOUT(IY) = VECT(IY) + (G1*VECT(IPY) + G2*HXP(2) + G3*F2)
         VOUT(IZ) = VECT(IZ) + (G1*VECT(IPZ) + G2*HXP(3) + G3*F3)
 
         VOUT(IPX) = VECT(IPX) + (G4*VECT(IPX) + G5*HXP(1) + G6*F1)
         VOUT(IPY) = VECT(IPY) + (G4*VECT(IPY) + G5*HXP(2) + G6*F2)
         VOUT(IPZ) = VECT(IPZ) + (G4*VECT(IPZ) + G5*HXP(3) + G6*F3)
*
      ELSE
         VOUT(IX) = VECT(IX) + STEP*VECT(IPX)
         VOUT(IY) = VECT(IY) + STEP*VECT(IPY)
         VOUT(IZ) = VECT(IZ) + STEP*VECT(IPZ)
*
      ENDIF
*
  999 END
