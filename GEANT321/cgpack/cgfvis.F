*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:43  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGFVIS(NT,FACE,IVIS,ISHAPE)
************************************************************************
*                                                                      *
*     Name: CGFVIS                                                     *
*     Authors: E. Chernyaev, S. Giani            Date:    01.08.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Transfer face to screen coordinates                    *
*               and find face min-max and face visibility              *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:  NT - number of transformation to screen coordinates      *
*             FACE - face                                              *
*                                                                      *
*     Output: IVIS - visibility flag                                   *
*                    1 - if visible face                               *
*                   -1 - if unvisible                                  *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cgcfac.inc"
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/cgctra.inc"
#include "geant321/gcspee.inc"
***SG
      DIMENSION ACCMI1(6),ACCMI2(6),ACCMI3(6),
     +          ACCMA1(6),ACCMA2(6),ACCMA3(6)
***SG
      REAL      FACE(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION  T(4,3),A,B,C,S
#endif
#if defined(CERNLIB_SINGLE)
      REAL      T(4,3)
#endif
*-
      IVIS   = -1
      DO 120 I=1,4
         DO 110 J=1,3
            T(I,J) = TSCRN(I,J,NT)
  110    CONTINUE
  120 CONTINUE
*
***SG
**         HIDDEN FACE REMOVAL
*     Computing face scope and skipping if it's 'covered': this
*     can allow a great increase in speed and a great reduction
*     in number of memory words used.
*
      IF(ISHAPE.LE.4.OR.ISHAPE.EQ.10)THEN
         NTIM=NTIM+1
         J = LCGFAC
         SRFMI1 = FACE(J+KCGX1)
         SRFMI2 = FACE(J+KCGY1)
         SRFMI3 = FACE(J+KCGZ1)
         SRFMA1 = FACE(J+KCGX1)
         SRFMA2 = FACE(J+KCGY1)
         SRFMA3 = FACE(J+KCGZ1)
         NEDGE = FACE(KCGNE)
         DO 333 NE=1,NEDGE
            SRFMI1 = MIN(SRFMI1,FACE(J+KCGX1),FACE(J+KCGX2))
            SRFMI2 = MIN(SRFMI2,FACE(J+KCGY1),FACE(J+KCGY2))
            SRFMI3 = MIN(SRFMI3,FACE(J+KCGZ1),FACE(J+KCGZ2))
            SRFMA1 = MAX(SRFMA1,FACE(J+KCGX1),FACE(J+KCGX2))
            SRFMA2 = MAX(SRFMA2,FACE(J+KCGY1),FACE(J+KCGY2))
            SRFMA3 = MAX(SRFMA3,FACE(J+KCGZ1),FACE(J+KCGZ2))
            J = J + LCGEDG
  333    CONTINUE
         IF(IPORLI.EQ.1)THEN
            ACCMI1(NTIM) = SRFMI1
            ACCMI2(NTIM) = SRFMI2
            ACCMI3(NTIM) = SRFMI3
            ACCMA1(NTIM) = SRFMA1
            ACCMA2(NTIM) = SRFMA2
            ACCMA3(NTIM) = SRFMA3
         ELSEIF(ISUBLI.EQ.1)THEN
            IF(NTIM.EQ.1.OR.NTIM.EQ.2)THEN
               IF(SRFMI3.LT.ACCMI3(NTIM).AND.SRFMA3
     +         .LT.ACCMA3(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.3.OR.NTIM.EQ.5)THEN
               IF(SRFMI2.LT.ACCMI2(NTIM).AND.SRFMA2
     +         .LT.ACCMA2(NTIM))GOTO 999
            ELSEIF(NTIM.EQ.4.OR.NTIM.EQ.6)THEN
               IF(SRFMI1.LT.ACCMI1(NTIM).AND.SRFMA1
     +         .LT.ACCMA1(NTIM))GOTO 999
            ENDIF
         ENDIF
      ENDIF
*
***SG
*
      C      = (T(2,1)*T(3,2) - T(3,1)*T(2,2))*FACE(KCGAA) +
     +         (T(3,1)*T(1,2) - T(1,1)*T(3,2))*FACE(KCGBB) +
     +         (T(1,1)*T(2,2) - T(2,1)*T(1,2))*FACE(KCGCC)
      IF (C .LE. 0.)     GOTO 999
      B      = (T(2,3)*T(3,1) - T(3,3)*T(2,1))*FACE(KCGAA) +
     +         (T(3,3)*T(1,1) - T(1,3)*T(3,1))*FACE(KCGBB) +
     +         (T(1,3)*T(2,1) - T(2,3)*T(1,1))*FACE(KCGCC)
      A      = (T(2,2)*T(3,3) - T(3,2)*T(2,3))*FACE(KCGAA) +
     +         (T(3,2)*T(1,3) - T(1,2)*T(3,3))*FACE(KCGBB) +
     +         (T(1,2)*T(2,3) - T(2,2)*T(1,3))*FACE(KCGCC)
      S      = 1./SQRT(A*A+B*B+C*C)
      AABCD(1) = A*S
      AABCD(2) = B*S
      AABCD(3) = C*S
*
      F1(KCGAF) = FACE(KCGAF)
      F1(KCGNE) = FACE(KCGNE)
      F1(KCGAA) = 0.
      F1(KCGBB) = 0.
      F1(KCGCC) = 1.
      F1(KCGDD) = 0.
      F1(KCGNE) = FACE(KCGNE)
*
**           T R A S F E R   P O I N T   C O O R D I N A T E S
*
      NEDGE  = FACE(KCGNE)
      IF (LCGFAC+NEDGE*LCGEDG .GT. LABC)
     +       PRINT *, ' Problem in CGFVIS: no space'
      XGRAV  = 0.
      YGRAV  = 0.
      ZGRAV  = 0.
      J      = LCGFAC
      DO 200 NE=1,NEDGE
         F1(J+KCGAE) = FACE(J+KCGAE)
         X = FACE(J+KCGX1)
         Y = FACE(J+KCGY1)
         Z = FACE(J+KCGZ1)
         F1(J+KCGX1) = T(1,1)*X + T(2,1)*Y + T(3,1)*Z + T(4,1)
         F1(J+KCGY1) = T(1,2)*X + T(2,2)*Y + T(3,2)*Z + T(4,2)
         F1(J+KCGZ1) = T(1,3)*X + T(2,3)*Y + T(3,3)*Z + T(4,3)
         X = FACE(J+KCGX2)
         Y = FACE(J+KCGY2)
         Z = FACE(J+KCGZ2)
         F1(J+KCGX2) = T(1,1)*X + T(2,1)*Y + T(3,1)*Z + T(4,1)
         F1(J+KCGY2) = T(1,2)*X + T(2,2)*Y + T(3,2)*Z + T(4,2)
         F1(J+KCGZ2) = T(1,3)*X + T(2,3)*Y + T(3,3)*Z + T(4,3)
         XGRAV = XGRAV + F1(J+KCGX1) + F1(J+KCGX2)
         YGRAV = YGRAV + F1(J+KCGY1) + F1(J+KCGY2)
         ZGRAV = ZGRAV + F1(J+KCGZ1) + F1(J+KCGZ2)
         J = J + LCGEDG
  200 CONTINUE
      XFACT = 1./(2.*NEDGE)
      XGRAV = XGRAV * XFACT
      YGRAV = YGRAV * XFACT
      ZGRAV = ZGRAV * XFACT
      AABCD(4) =-(AABCD(1)*XGRAV + AABCD(2)*YGRAV + AABCD(3)*ZGRAV)
*
**           F I N D   F A C E   M I N - M A X
*
      J      = LCGFAC
      RFMIN(1) = F1(J+KCGX1)
      RFMIN(2) = F1(J+KCGY1)
      RFMIN(3) = F1(J+KCGZ1)
      RFMAX(1) = F1(J+KCGX1)
      RFMAX(2) = F1(J+KCGY1)
      RFMAX(3) = F1(J+KCGZ1)
      DO 300 NE=1,NEDGE
         RFMIN(1) = MIN(RFMIN(1),F1(J+KCGX1),F1(J+KCGX2))
         RFMIN(2) = MIN(RFMIN(2),F1(J+KCGY1),F1(J+KCGY2))
         RFMIN(3) = MIN(RFMIN(3),F1(J+KCGZ1),F1(J+KCGZ2))
         RFMAX(1) = MAX(RFMAX(1),F1(J+KCGX1),F1(J+KCGX2))
         RFMAX(2) = MAX(RFMAX(2),F1(J+KCGY1),F1(J+KCGY2))
         RFMAX(3) = MAX(RFMAX(3),F1(J+KCGZ1),F1(J+KCGZ2))
         F1(J+KCGZ1) = 0.
         F1(J+KCGZ2) = 0.
         J = J + LCGEDG
  300 CONTINUE
      DRFACE(1) =-RFMAX(1)
      DRFACE(2) =-RFMAX(2)
      DRFACE(3) = RFMIN(1)
      DRFACE(4) = RFMIN(2)
      DRFACE(5) = RFMIN(3)
*
**          C O M P U T E    F A C E    V I S I B L E    A R E A
*
      J      = LCGFAC
      S      = 0.
      DLMAX  = 0.
      DO 400 NE=1,NEDGE
         S = S + F1(J+KCGX1)*F1(J+KCGY2) - F1(J+KCGX2)*F1(J+KCGY1)
         DL = ABS(F1(J+KCGX2)-F1(J+KCGX1)) + ABS(F1(J+KCGY2)-F1(J+
     +   KCGY1))
         IF (DLMAX .LT. DL) DLMAX = DL
         J = J + LCGEDG
  400 CONTINUE
      IF (DLMAX .LT. EESCR)     GOTO 999
      IF (S .GT. DLMAX*EESCR)   IVIS = 1
*
  999 RETURN
      END
