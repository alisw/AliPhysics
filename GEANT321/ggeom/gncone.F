*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:51  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GNCONE(X,P,IACT,IFL,SNEXT,SNXT,SAFE)
C.    ******************************************************************
C.    *                                                                *
C.    *      Compute distance to intersection with boundary surface of *
C     *      volume CONE or CONS, from point X(1),X(2),X(3) inside     *
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
C.    *      Authors   : Michel Maire and Rolf Nierhaus    21-JUN-1990 *
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
      REAL SNEXT, SNXT, SAFE
      PARAMETER (F=0.01745329251994330D0)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (F=0.01745329251994330)
#endif
      REAL X(6),P(7)
*
*     this part has to be moved outside the routine
      RO1=0.5*(P(4)+P(2))
      TG1=0.5*(P(4)-P(2))/P(1)
      CR1=1./SQRT(1.+TG1*TG1)
      RO2=0.5*(P(5)+P(3))
      TG2=0.5*(P(5)-P(3))/P(1)
      CR2=1./SQRT(1.+TG2*TG2)
      IF (IFL.EQ.2) THEN
         P6=P(6)*F
         P7=P(7)*F
         IF (P7.LT.P6) P7=P7+F*360.
         C1=COS(P6)
         S1=SIN(P6)
         C2=COS(P7)
         S2=SIN(P7)
         FIO=0.5*(P7+P6)
         CFIO=COS(FIO)
         SFIO=SIN(FIO)
      END IF
*
      SNXT=1.E10
      R   =SQRT(X(1)**2+X(2)**2)
      RIN =TG1*X(3)+RO1
      ROUT=TG2*X(3)+RO2
*
*     Compute SAFE radius
      IF (IACT.LT.3) THEN
         SAF1=(R -RIN)*CR1
         SAF2=(ROUT-R)*CR2
         SAF3=P(1)-ABS(X(3))
         SAF4=1.E10
         IF (IFL.EQ.2) THEN
            IF ((X(2)*CFIO-X(1)*SFIO).LE.0.) THEN
               SAF4=ABS(X(1)*S1-X(2)*C1)
            ELSE
               SAF4=ABS(X(1)*S2-X(2)*C2)
            END IF
         END IF
         SAFE=MIN(SAF1,SAF2,SAF3,SAF4)
         IF (IACT.EQ.0) GO TO 999
         IF (IACT.EQ.1.AND.SNEXT.LE.SAFE) GO TO 999
      END IF
*
*     Intersection with z-plane
      IF (X(6).GT. 1.E-20) THEN
         SZ= (P(1)-X(3))/X(6)
      ELSEIF (X(6).LT.-1.E-20) THEN
         SZ=-(P(1)+X(3))/X(6)
      ELSE
         SZ= 1.E10
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
*     Intersection with the inner cone
      SR1=1.E10
      IF (RO1.GT.0.) THEN
         U=T1-(TG1*X(6))**2
         V=T2- TG1*X(6)*RIN
         W=T3- RIN*RIN
*        track not parallel to the cone ?
         IF (U.NE.0.) THEN
            B=V/U
            C=W/U
            D=B**2-C
            IF (D.GE.0.) THEN
               DS=SQRT(D)
               IF (DS.GE.ABS(B)) THEN
                  SR1= DS-B
               ELSEIF (B.LE.0.)  THEN
                  SR1=-DS-B
               END IF
            END IF
         ELSEIF (V.LT.0.) THEN
            SR1=-0.5*W/V
         END IF
      END IF
*
*     Intersection with the outer cone
      SR2=1.E10
      U=T1-(TG2*X(6))**2
      V=T2- TG2*X(6)*ROUT
      W=T3- ROUT*ROUT
*     track not parallel to the cone ?
      IF (U.NE.0.) THEN
         B=V/U
         C=W/U
         D=B**2-C
         IF (D.GE.0.) THEN
            DS=SQRT(D)
            IF (DS.GE.ABS(B)) THEN
               SR2= DS-B
            ELSEIF (B.LE.0.)  THEN
               SR2=-DS-B
            END IF
         END IF
      ELSEIF (V.GT.0.) THEN
         SR2=-0.5*W/V
      END IF
*
*     Intersection with phi-planes
*     x=r*cos(phi)=X(1)+t*X(4)
*     y=r*sin(phi)=X(2)+t*X(5)
*     z           =X(3)+t*X(6)
*     t=(X(2)*cos(phi)-X(1)*sin(phi))/(X(4)*sin(phi)-X(5)*cos(phi))
*
      SFI1=1.E10
      SFI2=1.E10
      IF (IFL.EQ.2) THEN
*        track not parallel to the phi1 plane ?
         UN=X(4)*S1-X(5)*C1
         IF (UN.NE.0.) THEN
            S=(X(2)*C1-X(1)*S1)/UN
            IF (S.GE.0.) THEN
               XI=X(1)+S*X(4)
               YI=X(2)+S*X(5)
               IF ((YI*CFIO-XI*SFIO).LE.0.) SFI1=S
            END IF
         END IF
*        track not parallel to the phi2 plane ?
         UN=X(4)*S2-X(5)*C2
         IF (UN.NE.0.) THEN
            S=(X(2)*C2-X(1)*S2)/UN
            IF (S.GE.0.) THEN
               XI=X(1)+S*X(4)
               YI=X(2)+S*X(5)
               IF ((YI*CFIO-XI*SFIO).GE.0.) SFI2=S
            END IF
         END IF
      END IF
*
      SNXT=MIN(SR1,SR2,SZ,SFI1,SFI2)
  999 END
 
