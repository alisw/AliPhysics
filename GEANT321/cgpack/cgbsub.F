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
      SUBROUTINE CGBSUB(IFB,A,B,NMAX,C,IREP)
************************************************************************
*                                                                      *
*     Name: CGBSUB                                                     *
*     Author: E. Chernyaev                       Date:    23.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Subtraction of two faces                               *
*                                                                      *
*     References: CGBMMT,CGBTEF,CGBSOR                                 *
*                                                                      *
*     Input: IFB - flag (=1, if -B)
*            A(*) - 1-st face                                          *
*            B(*) - 2-nd face                                          *
*            NMAX - max length C array                                 *
*                                                                      *
*     Output: B(*) - resulting face                                    *
*             C(*) - working face                                      *
*             IREP - reply (legth of resulting face)                   *
*                          (-1 if no space)                            *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      A(*),B(*),C(*)
      DATA      IEXTNL/1/,IINTNL/-1/
*-
      NAEDGE  = A(KCGNE)
      NBEDGE  = B(KCGNE)
      IF (NBEDGE .EQ. 0)        GOTO 100
*           M I N - M A X   T E S T
      CALL CGBMMT(A,B,IREP)
      IF (IREP .NE. 0)          GOTO 200
*
**           N O   I N T E R S E C T I O N   O F   F A C E S
*
  100 IF (IFB .NE. 0)           GOTO 997
      LENB   = LCGFAC + NBEDGE*LCGEDG
      IREP   = LCGFAC + NAEDGE*LCGEDG
      IF (IREP .GT. NMAX+LENB)  GOTO 998
      DO 110 I=1,IREP
        B(I)   = A(I)
  110   CONTINUE
      GOTO 999
*
**          T E S T   1 - S T   F A C E   A G A I N S T   2 - N D
*
  200 IF (NMAX .LT. LCGFAC)     GOTO 998
      C(KCGAF) = NMAX
      C(KCGNE) = 0.
      IF (IFB .EQ. 0)           IFLAG = IEXTNL
      IF (IFB .NE. 0)           IFLAG = IINTNL
      IVAR   = 1
      J      = LCGFAC
      DO 210 NE=1,NAEDGE
        CALL CGBTEF(IFLAG,IVAR,A(J+1),B,C)
        IF (C(KCGAF) .LT. 0.)   GOTO 998
        J      = J + LCGEDG
  210   CONTINUE
*
**           T E S T   2 - N D   F A C E   A G A I N S T   1 - S T
*
      J      = LCGFAC
      DO 300 NE=1,NBEDGE
        CALL CGBTEF(IINTNL,IVAR,B(J+1),A,C)
        IF (C(KCGAF) .LT. 0.)   GOTO 998
        J      = J + LCGEDG
  300   CONTINUE
*           P R E P A R E   N E W   F A C E
      NEDGE = C(KCGNE)
      IF (NEDGE .EQ. 0)         GOTO 997
      CALL CGBSOR(NEDGE,C(LCGFAC+1))
      C(KCGNE) = NEDGE
      IF (NEDGE .EQ. 0)         GOTO 997
*      IF (NEDGE .LT. 3) PRINT *,' CGBSUB: NEDGE .LT. 3 - face ignored'
      C(KCGAF) = A(KCGAF)
      C(KCGAA) = A(KCGAA)
      C(KCGBB) = A(KCGBB)
      C(KCGCC) = A(KCGCC)
      C(KCGDD) = A(KCGDD)
      IREP  = LCGFAC + NEDGE*LCGEDG
      DO 400 I=1,IREP
        B(I)   = C(I)
  400   CONTINUE
      GOTO 999
*
  997 IREP   = 0
      GOTO 999
  998 IREP = -1
  999 RETURN
      END
 
