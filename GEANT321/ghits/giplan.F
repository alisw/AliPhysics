*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:10  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.20  by  S.Giani
*-- Author :
      SUBROUTINE GIPLAN(YC,X1,X2,S1,S2,IC,XINT,SINT,PZINT,IFLAG)
C.
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Calculates intersection of track (X1,X2)                 *
C.    *       with plane parallel to (X-Z)                             *
C.    *        The track is approximated by a cubic in the             *
C.    *       track length.                                            *
C.    *       To improve stability, the coordinate system              *
C.    *       is shifted.                                              *
C.    *       input parameters                                         *
C.    *        YC    = Y COORDINATE OF PLANE                           *
C.    *        X1    = X,Y,Z,XP,YP,ZP OF 1ST POINT                     *
C.    *        X2    =                   2ND                           *
C.    *        S1(2) = S AT 1ST(2ND) POINT                             *
C.    *        IC    = 1 STRAIGHT LINE DEFINED BY X+XP                 *
C.    *        IC    = 2 STRAIGHT LINE DEFINED BY X1+X2                *
C.    *        IC    = 3 CUBIC MODEL                                   *
C.    *                                                                *
C.    *      output parameters                                         *
C.    *        XINT  = X,Y,Z,XP,YP,ZP AT INTERSECTION POINT            *
C.    *        SINT  = S AT INTERSECTION POINT                         *
C.    *        PZINT = PHI,Z,DPHI/DR,DZ/DR                             *
C.    *        IFLAG = 1 IF TRACK INTERSECTS PLANE                     *
C.    *              = 0 IF NOT                                        *
C.    *                                                                *
C.    *      Warning : the default accuracy is 10 microns. The value   *
C.    *      of EPSI must be changed for a better precision            *
C.    *                                                                *
C.    *    ==>Called by : <USER>, GUDIGI                               *
C.    *                                                                *
C.    *        Authors: R.BRUN/JJ.DUMONT from an original routine by   *
C.    *       H. BOERNER  KEK  OCTOBER 1982                            *
C.    *                                                                *
C.    *                                                                *
C.    ******************************************************************
C.
      DIMENSION X1(6),X2(6),XINT(6),PZINT(4),A(4),B(4),C(4)
C
      DATA MAXHIT/100/
      DATA EPSI/0.001/
C.
C.
C.    ------------------------------------------------------------------
C.
C.
      IFLAG  = 1
      DRCTN  = 1.
C
C             Track crossing the plane from above or below ?
C
      IF (X2(2).LT.X1(2))                        GO TO 5
      IF (YC.LT.X1(2))                           GO TO 90
      IF (YC.GT.X2(2))                           GO TO 90
      IF (IC.EQ.2) GO TO 30
      IF (IC.EQ.3) GO TO 7
C
      S=S1
      DXDS=X1(4)
      DYDS=X1(5)
      DZDS=X1(6)
      BX=X1(1)-DXDS*(X1(2)-YC)/DYDS
      BZ=X1(3)-DZDS*(X1(2)-YC)/DYDS
      TRL2=(BX-X1(1))**2+(X1(2)-YC)**2+(BZ-X1(3))**2
      GO TO 40
C
   5  IF (YC.LT.X2(2))                           GO TO 90
      IF (YC.GT.X1(2))                           GO TO 90
      IF(IC.EQ.2) GO TO 30
      DRCTN  = - 1.
C
      IF(IC.EQ.3) GOTO 7
      S=S2
      DXDS=X2(4)
      DYDS=X2(5)
      DZDS=X2(6)
      BX=X2(1)-DXDS*(X2(2)-YC)/DYDS
      BZ=X2(3)-DZDS*(X2(2)-YC)/DYDS
      TRL2=(BX-X2(1))**2+(X2(2)-YC)**2+(BZ-X2(3))**2
      GOTO 40
C
   30 DX=X2(1)-X1(1)
      DY=X2(2)-X1(2)
      DZ=X2(3)-X1(3)
      DS=SQRT(DX*DX+DY*DY+DZ*DZ)
      S=S1
      DXDS=DX/DS
      DYDS=DY/DS
      DZDS=DZ/DS
      BX=X1(1)-DX*(X1(2)-YC)/DY
      BZ=X1(3)-DZ*(X1(2)-YC)/DY
      TRL2=(BX-X1(1))**2+(X1(2)-YC)**2+(BZ-X1(3))**2
C
   40 TRLEN=SQRT(TRL2)*DRCTN+S
      XINT(1)=BX
      XINT(2)=YC
      XINT(3)=BZ
      SINT=TRLEN
      XINT(4)=DXDS
      XINT(5)=DYDS
      XINT(6)=DZDS
      GO TO 200
C
C               Shift coordinate system such that center of gravity=0
C
   7  IF(YC.LE.0.)                                GO TO 90
      SHIFTX = (X1(1) + X2(1)) * 0.5
      SHIFTY = (X1(2) + X2(2)) * 0.5
      SHIFTZ = (X1(3) + X2(3)) * 0.5
      SHIFTS = (S1 + S2) * 0.5
C
C             Only one value necessary since X1= -X2 etc...
C
      XSHFT  = X1(1) - SHIFTX
      YSHFT  = X1(2) - SHIFTY
      ZSHFT  = X1(3) - SHIFTZ
      SSHFT  = S1 - SHIFTS
C
      PABS1  = SQRT(X1(4)**2 + X1(5)**2 + X1(6)**2)
      PABS2  = SQRT(X2(4)**2 + X2(5)**2 + X2(6)**2)
      IF (PABS1.EQ.0..OR.PABS2.EQ.0.)            GO TO 90
C
C              Parametrize the track by a cubic through X1, X2
C
      CALL GCUBS(SSHFT,XSHFT,X1(4)/PABS1,X2(4)/PABS2,A)
      CALL GCUBS(SSHFT,YSHFT,X1(5)/PABS1,X2(5)/PABS2,B)
      CALL GCUBS(SSHFT,ZSHFT,X1(6)/PABS1,X2(6)/PABS2,C)
C
C              Iterate to find the track length corresponding to
C              the intersection of track and plane.
C              Start at S=0. middle of the shifted interval.
C
      DINTER = ABS(S2 - S1) * 0.5
      S      = 0.
C
      DO 10 I = 1,MAXHIT
         Y = SHIFTY + B(1) + S * (B(2) + S * (B(3) + S * B(4)))
         DR=(YC-Y)*DRCTN
         IF (ABS(DR).LT.EPSI)                     GO TO 20
         DINTER = DINTER * 0.5
         IF (DR.LT.0.)S = S - DINTER
         IF (DR.GE.0.)S = S + DINTER
  10  CONTINUE
C
C              Compute intersection in original coordinates
C
  20  CONTINUE
      XINT(1) = SHIFTX + A(1) + S * (A(2) + S * (A(3) + S * A(4)))
      XINT(2)=YC
      XINT(3) = SHIFTZ + C(1) + S * (C(2) + S * (C(3) + S * C(4)))
      XINT(4) = A(2) + S * (2. * A(3) + 3. * S * A(4))
      XINT(5) = B(2) + S * (2. * B(3) + 3. * S * B(4))
      XINT(6) = C(2) + S * (2. * C(3) + 3. * S * C(4))
C
C             Compute PHIHIT,ZHIT and corresponding derivatives
C
      SINT   = S + SHIFTS
  200 TERM   = 1. / (XINT(4) * XINT(1) + XINT(5) * XINT(2))
      PZINT(1) = ATAN2(XINT(2),XINT(1))
      PZINT(2) = XINT(3)
      PZINT(3) = (XINT(1) * XINT(5) - XINT(2) * XINT(4)) * TERM / YC
      PZINT(4) = TERM * XINT(6) * YC
      RETURN
C
  90  IFLAG  = 0
      END
