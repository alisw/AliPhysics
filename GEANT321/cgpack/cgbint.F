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
      SUBROUTINE CGBINT(A,B,NMAX,C,IREP)
************************************************************************
*                                                                      *
*     Name: CGBINT                                                     *
*     Author: E. Chernyaev                       Date:    23.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Intersection of two faces                              *
*                                                                      *
*     References: CGBMMT,CGBTEF,CGBSOR                                 *
*                                                                      *
*     Input: A(*) - 1-st face                                          *
*            B(*) - 2-nd face                                          *
*            NMAX - max length C array                                 *
*                                                                      *
*     Output: C(*) - resulting face                                    *
*             IREP - reply (legth of resulting face)                   *
*                          (-1 if no space)                            *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      A(*),B(*),C(*)
*-
      IF (NMAX .LT. LCGFAC)     GOTO 998
*           M I N - M A X   T E S T   O F   T W O   F A C E S
      CALL CGBMMT(A,B,IREP)
      IF (IREP .EQ. 0)          GOTO 999
*           T E S T   1 - S T   F A C E   A G A I N S T   2 - N D
      C(KCGAF) = NMAX
      C(KCGNE) = 0.
      NEDGE  = A(KCGNE)
      J      = LCGFAC
      INTRNL = -1
      IVAR   = 2
      DO 100 NE=1,NEDGE
        CALL CGBTEF(INTRNL,IVAR,A(J+1),B,C)
        IF (C(KCGAF) .LT. 0.)   GOTO 998
        J      = J + LCGEDG
  100   CONTINUE
*           T E S T   2 - N D   F A C E   A G A I N S T   1 - S T
      NEDGE  = B(KCGNE)
      J      = LCGFAC
      DO 200 NE=1,NEDGE
        CALL CGBTEF(INTRNL,IVAR,B(J+1),A,C)
        IF (C(KCGAF) .LT. 0.)   GOTO 998
        J      = J + LCGEDG
  200   CONTINUE
*           P R E P A R E   N E W   F A C E
      NEDGE = C(KCGNE)
      IF (NEDGE .EQ. 0)         GOTO 997
      CALL CGBSOR(NEDGE,C(LCGFAC+1))
      C(KCGNE) = NEDGE
      IF (NEDGE .EQ. 0)         GOTO 997
      IF (NEDGE .LT. 3) PRINT *,' CGBINT: NEDGE .LT. 3 - face ignored'
      IF (NEDGE .LT. 3)         GOTO 997
      C(KCGAF) = A(KCGAF)
      IF (B(KCGAF) .GT. A(KCGAF)) C(KCGAF) = B(KCGAF)
      C(KCGAA) = A(KCGAA)
      C(KCGBB) = A(KCGBB)
      C(KCGCC) = A(KCGCC)
      C(KCGDD) = A(KCGDD)
      IREP  = LCGFAC + NEDGE*LCGEDG
      GOTO 999
*
  997 IREP   = 0
      GOTO 999
  998 IREP = -1
  999 RETURN
      END
