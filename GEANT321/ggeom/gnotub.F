*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:53  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.30  by  S.Giani
*-- Author :
      SUBROUTINE GNOTUB(X,P,IACT,IFL,SNEXT,SNXT,SAFE)
C.    ******************************************************************
C.    *                                                                *
C.    *      Compute distance to intersection with boundary surface of *
C     *      volume TUBE or TUBS, from point X(1),X(2),X(3) outside    *
C     *      the volume along track with direction cosines X(4),X(5),  *
C     *      X(6)                                                      *
C.    *      P     (input)  : volume parameters                        *
C.    *      IACT  (input)  : action flag                              *
C.    *         = 0  Compute SAFE only                                 *
C.    *         = 1  Compute SAFE, compute SNXT only if SAFE.LT.SNEXT  *
C.    *         = 2  Compute both SAFE and SNXT                        *
C.    *         = 3  Compute SNXT only                                 *
C.    *      IFL   (input)  : 1 for TUBE, 2 for PHI segmented TUBE     *
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
C.    * 'TUBE'    is a tube. It has 3 parameters, the inside radius,   *
C.    *           the outside radius and the half length in z.         *
C.    * 'TUBS'    is a phi segment of a  tube.  It has 5 parameters,   *
C.    *           the same  3 as  'TUBE' plus  the phi  limits.  The   *
C.    *           segment  starts at  the first  limit and  includes   *
C.    *           increasing phi  value up  to the  second limit  or   *
C.    *           that plus 360 degrees.                               *
C.    *                                                                *
C.    ******************************************************************
#if !defined(CERNLIB_SINGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#include "geant321/gconsp.inc"
      REAL X(6),P(5),SNEXT,SNXT,SAFE
*
*     this part has to be moved outside the routine
      IF (IFL.EQ.2) THEN
         P4=P(4)*DEGRAD
         P5=P(5)*DEGRAD
         IF (P5.LT.P4) P5=P5+TWOPI
         C1=COS(P4)
         S1=SIN(P4)
         C2=COS(P5)
         S2=SIN(P5)
         FIO=0.5*(P5+P4)
         CFIO=COS(FIO)
         SFIO=SIN(FIO)
         DFI=0.5*(P5-P4)
         CDFI=COS(DFI)
*        SDFI=SIN(DFI)
      END IF
*
      SNXT=1.E10
      R=SQRT(X(1)**2+X(2)**2)
*
*     Compute SAFE radius
      IF (IACT.LT.3) THEN
         SAF1=P(1)-R
         SAF2=R-P(2)
         SAF3=ABS(X(3))-P(3)
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
         ENDIF
         SAFE=MAX(SAF1,SAF2,SAF3,SAF4)
         IF (IACT.EQ.0) GO TO 999
         IF (IACT.EQ.1.AND.SNEXT.LE.SAFE) GO TO 999
      END IF
*
*     Intersection with z-plane
*     only points outside the z range need to be considered
      IF (ABS(X(3)).GE.P(3)) THEN
         IF (X(3)*X(6).LT.0.) THEN
            S=(ABS(X(3))-P(3))/ABS(X(6))
            XI=X(1)+S*X(4)
            YI=X(2)+S*X(5)
            R2=XI**2+YI**2
            IF (P(1)**2.LE.R2.AND.R2.LE.P(2)**2) THEN
               IF (IFL.EQ.1) GO TO 101
               IF (R2.EQ.0.) GO TO 101
               CPSI=(XI*CFIO+YI*SFIO)/SQRT(R2)
               IF (CPSI.GE.CDFI) GO TO 101
            END IF
         END IF
      END IF
*
*     Intersection with cylinders
*     Intersection point (x,y,z)
*     (x,y,z) is on track :    x=X(1)+t*X(4)
*                              y=X(2)+t*X(5)
*                              z=X(3)+t*X(6)
*     (x,y,z) is on cylinder : x**2 + y**2 = R**2
*
*     (X(4)**2+X(5)**2)*t**2
*     +2.*(X(1)*X(4)+X(2)*X(5))*t
*     +X(1)**2+X(2)**2-R**2=0
*
      T1=X(4)**2+X(5)**2
      T2=(X(1)*X(4)+X(2)*X(5))
      T3=X(1)**2+X(2)**2
*     track parallel to the z axis ?
***** 21-JUN-1990
*     IF (T1.EQ.0.) GO TO 999
      IF (ABS(T1).LT.1.E-32) GO TO 999
      B =T2/T1
*
*     Intersection with outer cylinder
*     only points outside the outer cylinder need to be considered
      IF (R.GE.P(2)) THEN
         C=(T3-P(2)**2)/T1
         D=B**2-C
         IF (D.GE.0.) THEN
            S=-B-SQRT(D)
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(3)) THEN
                  IF (IFL.EQ.1) GO TO 101
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  CPSI=(XI*CFIO+YI*SFIO)/P(2)
                  IF (CPSI.GE.CDFI) GO TO 101
               END IF
            END IF
         END IF
      END IF
*     Intersection with inner cylinder
      IF (P(1).GT.0.) THEN
         C=(T3-P(1)**2)/T1
         D=B**2-C
         IF (D.GE.0.) THEN
            S=-B+SQRT(D)
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(3)) THEN
                  IF (IFL.EQ.1) GO TO 101
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  CPSI=(XI*CFIO+YI*SFIO)/P(1)
                  IF (CPSI.GE.CDFI) SNXT=S
               END IF
            END IF
         END IF
      END IF
*
*     Intersection with phi-planes
*     x=r*cos(phi)=X(1)+t*X(4)
*     y=r*sin(phi)=X(2)+t*X(5)
*     z           =X(3)+t*X(6)
*     t=(X(2)*cos(phi)-X(1)*sin(phi))/(X(4)*sin(phi)-X(5)*cos(phi))
      IF (IFL.EQ.2) THEN
*        track not parallel to the phi1 plane ?
         UN=X(4)*S1-X(5)*C1
         IF (UN.NE.0.) THEN
            S=(X(2)*C1-X(1)*S1)/UN
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(3)) THEN
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  R2=XI**2+YI**2
                  IF (P(1)**2.LE.R2.AND.R2.LE.P(2)**2) THEN
                     IF ((YI*CFIO-XI*SFIO).LE.0.) THEN
                        IF (S.LT.SNXT) SNXT=S
                     END IF
                  END IF
               END IF
            END IF
         END IF
*        track not parallel to the phi2 plane ?
         UN=X(4)*S2-X(5)*C2
         IF (UN.NE.0.) THEN
            S=(X(2)*C2-X(1)*S2)/UN
            IF (S.GE.0.) THEN
               ZI=X(3)+S*X(6)
               IF (ABS(ZI).LE.P(3)) THEN
                  XI=X(1)+S*X(4)
                  YI=X(2)+S*X(5)
                  R2=XI**2+YI**2
                  IF (P(1)**2.LE.R2.AND.R2.LE.P(2)**2) THEN
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
 
