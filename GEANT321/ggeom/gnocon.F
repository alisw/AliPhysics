*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:52  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNOCON(X,P,IACT,IFL,SNEXT,SNXT,SAFE)
C.    ******************************************************************
C.    *                                                                *
C.    *      Compute distance to intersection with boundary surface of *
C     *      volume CONE or CONS, from point X(1),X(2),X(3) outside    *
C     *      the volume along track with direction cosines X(4),X(5),  *
C     *      X(6)                                                      *
C.    *      P     (input)  : volume parameters                        *
C.    *      IACT  (input)  : action flag                              *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, compute SNXT only if SAFE.LT.SNEXT  *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *      IFL   (input)  : 1 for CONE, 2 for PHI segmented CONE     *
C.    *      SNEXT (input)  : see IACT = 1                             *
C.    *      SNXT  (output) : distance to volume boundary along track  *
C.    *      SAFE  (output) : not larger than scalar distance to       *
C.    *                       volume boundaray                         *
C.    *      Called by : GNEXT, GNOPCO, GTNEXT                         *
C.    *                                                                *
C.    *      Authors   : Michel Maire and Rolf Nierhaus   21-JUN-1990  *
C.    *                                                                *
C.    ******************************************************************
C.    *                                                                *
C.    * 'CONE'    is a conical tube. It has 5 parameters :             *
C.    *           the half length in z,                                *
C.    *           the inside and outside radii at the low z limit,     *
C.    *           and those at the high z limit.                       *
C.    * 'CONS'    is a phi segment of a  conical tube.  It has 7       *
C.    *           parameters, the same 5 as 'CONE' plus the phi limits.*
C.    *           The segment starts at the first limit and  includes  *
C.    *           increasing phi  value up  to the  second limit  or   *
C.    *           that plus 360 degrees.                               *
C.    *                                                                *
C.    ******************************************************************
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (F=0.01745329251994330D0)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (F=0.01745329251994330)
#endif
      REAL X(6),P(7),SNEXT,SNXT,SAFE
      PARAMETER (ONE=1,HALF=ONE/2,ZERO=0)
*
*     this part has to be moved outside the routine
      RO1=HALF*(P(4)+P(2))
      TG1=HALF*(P(4)-P(2))/P(1)
      CR1=ONE/SQRT(ONE+TG1*TG1)
      ZV1=1.E10
      IF (P(2).NE.P(4)) ZV1=-RO1/TG1
      RO2=HALF*(P(5)+P(3))
      TG2=HALF*(P(5)-P(3))/P(1)
      CR2=ONE/SQRT(ONE+TG2*TG2)
      ZV2=1.E10
      IF (P(3).NE.P(5)) ZV2=-RO2/TG2
      IF (IFL.EQ.2) THEN
         P6=P(6)*F
         P7=P(7)*F
         IF (P7.LT.P6) P7=P7+F*360
         C1=COS(P6)
         S1=SIN(P6)
         C2=COS(P7)
         S2=SIN(P7)
         FIO=HALF*(P7+P6)
         CFIO=COS(FIO)
         SFIO=SIN(FIO)
         DFI=HALF*(P7-P6)
         CDFI=COS(DFI)
*        SDFI=SIN(DFI)
      END IF
*
      SNXT=1.E10
      R   =SQRT(X(1)**2+X(2)**2)
      RIN =ABS(TG1*X(3)+RO1)
      ROUT=ABS(TG2*X(3)+RO2)
*
*     Compute SAFE radius
      IF (IACT.LT.3) THEN
         SAF1=(RIN -R)*CR1
         SAF2=(R-ROUT)*CR2
         SAF3=ABS(X(3))-P(1)
         SAF4=0.
         IF (IFL.EQ.2.AND.R.GT.0.) THEN
            CPSI=(X(1)*CFIO+X(2)*SFIO)/R
            IF (CPSI.LT.CDFI) THEN
               IF ((X(2)*CFIO-X(1)*SFIO).LE.0.) THEN
                  SAF4=ABS(X(1)*S1-X(2)*C1)
               ELSE
                  SAF4=ABS(X(1)*S2-X(2)*C2)
               END IF
            END IF
         END IF
         SAFE=MAX(SAF1,SAF2,SAF3,SAF4)
         IF (IACT.EQ.0) GO TO 999
         IF (IACT.EQ.1.AND.SNEXT.LE.SAFE) GO TO 999
      END IF
*
*     Intersection with z-plane
*     only points outside the z range need to be considered
      IF (ABS(X(3)).GE.P(1)) THEN
         IF (X(3)*X(6).LT.0.) THEN
            S=(ABS(X(3))-P(1))/ABS(X(6))
            XI=X(1)+S*X(4)
            YI=X(2)+S*X(5)
            RIQ=XI**2+YI**2
            IF (X(3).LE.0.) THEN
               R1Q=P(2)**2
               R2Q=P(3)**2
            ELSE
               R1Q=P(4)**2
               R2Q=P(5)**2
            END IF
            IF (R1Q.LE.RIQ.AND.RIQ.LE.R2Q) THEN
               IF (IFL.EQ.1.OR.RIQ.LE.0.)  GO TO 101
               CPSI=(XI*CFIO+YI*SFIO)/SQRT(RIQ)
               IF (CPSI.GE.CDFI) GO TO 101
            END IF
         END IF
      END IF
*
*     Intersection with cones
*     Intersection point (x,y,z)
*     (x,y,z) is on track: x=X(1)+t*X(4)
*                          y=X(2)+t*X(5)
*                          z=X(3)+t*X(6)
*     (x,y,z) is on cone : x**2 + y**2 = (a*z+b)**2
*
*     (X(4)**2+X(5)**2-(a*x(6))**2)*t**2
*     +2.*(X(1)*X(4)+X(2)*X(5)-a*x(6)*(a*x(3)+b))*t
*     +X(1)**2+X(2)**2-(a*x(3)+b)**2=0
*
      T1=X(4)**2+X(5)**2
      T2=X(1)*X(4)+X(2)*X(5)
      T3=X(1)**2+X(2)**2
*
*     Intersection with the outer cone
*     only points outside the outer cone need to be considered
      SR2=1.E10
      IF ((ZV2*X(3).GT.ZV2*ZV2).OR.(R.GT.ROUT)) THEN
         U=T1-(TG2*X(6))**2
         V=T2- TG2*X(6)*(TG2*X(3)+RO2)
         W=T3- ROUT*ROUT
*        track not parallel to the cone ?
         IF (U.NE.0.) THEN
            B=V/U
            C=W/U
            D=B**2-C
            IF (D.GE.0.) THEN
*              compute the smallest root first
               IR=-1
   11          S=-B+IR*SQRT(D)
               IF (S.GE.0.) THEN
                  ZI=X(3)+S*X(6)
                  IF (ABS(ZI).LE.P(1)) THEN
                     IF (IFL.EQ.1.OR.ZI.EQ.ZV2) THEN
                        SR2=S
                     ELSE
                        XI=X(1)+S*X(4)
                        YI=X(2)+S*X(5)
                        RI=TG2*ZI+RO2
                        CPSI=(XI*CFIO+YI*SFIO)/RI
                        IF (CPSI.GE.CDFI) SR2=S
                     END IF
                  END IF
               END IF
               IF (IR.EQ.-1) THEN
                  IF (SR2.EQ.S) GO TO 101
*                 smallest root not ok. Try the biggest one
                  IR=+1
                  GO TO 11
               END IF
            END IF
         ELSEIF (V.NE.0.) THEN
            S=-HALF*W/V
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(1)) THEN
                  IF (IFL.EQ.1.OR.ZI.EQ.ZV2)  GO TO 101
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  RI=TG2*ZI+RO2
                  CPSI=(XI*CFIO+YI*SFIO)/RI
                  IF (CPSI.GE.CDFI) GO TO 101
               END IF
            END IF
         END IF
      END IF
*
*     Intersection with the inner cone
      SR1=1.E10
      IF (RO1.GT.0.) THEN
         U=T1-(TG1*X(6))**2
         V=T2- TG1*X(6)*(TG1*X(3)+RO1)
         W=T3- RIN*RIN
*        track not parallel to the cone ?
         IF (U.NE.0.) THEN
            B=V/U
            C=W/U
            D=B**2-C
            IF (D.GE.0.) THEN
*              compute the smallest root first
               IR=-1
   21          S=-B+IR*SQRT(D)
               IF (S.GE.0.) THEN
                  ZI=X(3)+S*X(6)
                  IF (ABS(ZI).LE.P(1)) THEN
                     IF (IFL.EQ.1.OR.ZI.EQ.ZV1) THEN
                        SR1=S
                     ELSE
                        XI=X(1)+S*X(4)
                        YI=X(2)+S*X(5)
                        RI=TG1*ZI+RO1
                        CPSI=(XI*CFIO+YI*SFIO)/RI
                        IF (CPSI.GE.CDFI) SR1=S
                     END IF
                  END IF
               END IF
               IF (IR.EQ.-1.AND.SR1.GT.9.E9) THEN
*                 smallest root not ok. Try the biggest one
                  IR=+1
                  GO TO 21
               END IF
            END IF
         ELSEIF (V.NE.0.) THEN
            S=-HALF*W/V
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(1)) THEN
                  IF (IFL.EQ.1.OR.ZI.EQ.ZV1) THEN
                     SR1=S
                  ELSE
                     XI=X(1)+S*X(4)
                     YI=X(2)+S*X(5)
                     RI=TG1*ZI+RO1
                     CPSI=(XI*CFIO+YI*SFIO)/RI
                     IF (CPSI.GE.CDFI) SR1=S
                  END IF
               END IF
            END IF
         END IF
      END IF
*
      SNXT=MIN(SR1,SR2)
*
*     Intersection with phi-planes
*     x=r*cos(phi)=X(1)+t*X(4)
*     y=r*sin(phi)=X(2)+t*X(5)
*     z           =X(3)+t*X(6)
*     t=(X(2)*cos(phi)-X(1)*sin(phi))/(X(4)*sin(phi)-X(5)*cos(phi))
*
      IF (IFL.EQ.2) THEN
*        track not parallel to the phi1 plane ?
         UN=X(4)*S1-X(5)*C1
         IF (UN.NE.0.) THEN
            S=(X(2)*C1-X(1)*S1)/UN
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(1)) THEN
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  RIQ=XI**2+YI**2
                  R1Q=(TG1*ZI+RO1)**2
                  R2Q=(TG2*ZI+RO2)**2
                  IF (R1Q.LE.RIQ.AND.RIQ.LE.R2Q) THEN
                     IF ((YI*CFIO-XI*SFIO).LE.0.) THEN
                        IF (S.LT.SNXT) SNXT=S
                     END IF
                  END IF
 
              END IF
            END IF
         END IF
*        track not parallel to the phi2 plane ?
         UN=X(4)*S2-X(5)*C2
         IF(UN.NE.0.) THEN
            S=(X(2)*C2-X(1)*S2)/UN
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(1)) THEN
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  RIQ=XI**2+YI**2
                  R1Q=(TG1*ZI+RO1)**2
                  R2Q=(TG2*ZI+RO2)**2
                  IF (R1Q.LE.RIQ.AND.RIQ.LE.R2Q) THEN
                     IF ((YI*CFIO-XI*SFIO).GE.0.) THEN
                        IF (S.LT.SNXT) SNXT=S
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      GO TO 999
*
  101 SNXT=S
  999 END
