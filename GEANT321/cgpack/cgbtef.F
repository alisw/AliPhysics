*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:42  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGBTEF(IFEXT,IVAR,EDGE,FACE,C)
************************************************************************
*                                                                      *
*     Name: CGBTEF                                                     *
*     Author: E. Chernyaev                       Date:    23.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Test edge against face                                 *
*                                                                      *
*     References: CGBEDG, CGTSOR                                       *
*                                                                      *
*     Input: IFEXT - flag for searching external or internal intervals *
*                    ( 1 - external, -1 - internal)                    *
*             IVAR - number of variant                                 *
*                    ( 1 - visibility of intervals the same as edge)   *
*                    ( 2 - intervals are invisible)                    *
*          EDGE(*) - number of edges in face                           *
*        FACE(*,*) - face                                              *
*             NMAX - max length of C array                             *
*                                                                      *
*     Output: C(*) - a set of new edges                                *
*                    (C(KCGAF) = -1 if no space)                       *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgcedg.inc"
      CHARACTER*2 WHAT
      REAL      EDGE(LCGEDG),FACE(*),C(*),ABCD(4)
      INTEGER   KCG(6)
      DATA      KCG/KCGX1,KCGY1,KCGZ1,KCGX2,KCGY2,KCGZ2/
*-
**          F I N D   I N T E R S E C T I O N   P O I N T S
*
      IF (IFEXT .LE. 0)         WHAT = 'LE'
      IF (IFEXT .GT. 0)         WHAT = 'GE'
      XD     = EDGE(KCGX2) - EDGE(KCGX1)
      YD     = EDGE(KCGY2) - EDGE(KCGY1)
      ZD     = EDGE(KCGZ2) - EDGE(KCGZ1)
      ALENG  = SQRT(XD*XD + YD*YD + ZD*ZD)
      IF(ALENG.LT.1.0E-4)GOTO 998
      XD     = XD / ALENG
      YD     = YD / ALENG
      ZD     = ZD / ALENG
      ABCD(1)= YD*FACE(KCGCC) - FACE(KCGBB)*ZD
      ABCD(2)= ZD*FACE(KCGAA) - FACE(KCGCC)*XD
      ABCD(3)= XD*FACE(KCGBB) - FACE(KCGAA)*YD
      ABCD(4)=-(ABCD(1)*EDGE(KCGX1) +
     +          ABCD(2)*EDGE(KCGY1) +
     +          ABCD(3)*EDGE(KCGZ1))
      CALL CGBFIT(FACE,ABCD,NT)
      IF (NT .GT. 0)                    GOTO 100
      XA     = EDGE(KCGX1)
      YA     = EDGE(KCGY1)
      ZA     = EDGE(KCGZ1)
      XDELT  = EDGE(KCGX2) - EDGE(KCGX1)
      YDELT  = EDGE(KCGY2) - EDGE(KCGY1)
      ZDELT  = EDGE(KCGZ2) - EDGE(KCGZ1)
*
**          P R E P A R E   E D G E S
*
  100 K      = 1
      IF (ABS(YDELT) .GT. ABS(XDELT))   K = 2
      IF (ABS(ZDELT) .GT. ABS(DELTA(K))) K = 3
      TMIN   = (EDGE(KCG(K))    - AA(K)) / DELTA(K)
      TMAX   = (EDGE(KCG(K+3))  - AA(K)) / DELTA(K)
      CALL CGBTTT(WHAT,TMIN,TMAX,NT,NEDGE)
      IF (NEDGE .EQ. 0)                 GOTO 999
      IVIS   = EDGE(KCGAE)
      NMAX   = C(KCGAF)
      NN     = C(KCGNE)
      J      = LCGFAC + NN*LCGEDG
      IF (NMAX .LT. J + NEDGE*LCGEDG)   GOTO 998
      DO 500 NE=1,NEDGE
        IF (IVAR .EQ. 1)                C(J+KCGAE) = IVIS
        IF (IVAR .EQ. 2)                C(J+KCGAE) =-1.
        IF (ITTT(NE) .NE. 0)            C(J+KCGAE) = IVIS
        C(J+KCGX1) = XA + XDELT*TTT(1,NE)
        C(J+KCGY1) = YA + YDELT*TTT(1,NE)
        C(J+KCGZ1) = ZA + ZDELT*TTT(1,NE)
        C(J+KCGX2) = XA + XDELT*TTT(2,NE)
        C(J+KCGY2) = YA + YDELT*TTT(2,NE)
        C(J+KCGZ2) = ZA + ZDELT*TTT(2,NE)
        J      = J + LCGEDG
  500   CONTINUE
      C(KCGNE) = C(KCGNE) + NEDGE
      GOTO 999
*
  998 C(KCGAF) =-1.
  999 RETURN
      END
