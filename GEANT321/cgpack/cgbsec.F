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
      SUBROUTINE CGBSEC(FACE,CG,RMN,RMX,NMAX,C,IREP)
************************************************************************
*                                                                      *
*     Name: CGBTOP                                                     *
*     Author: E. Chernyaev                       Date:    24.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Section of CG-object by face plane                     *
*                                                                      *
*     References: CGBTFP, CGBSOR                                       *
*                                                                      *
*     Input: FACE(*) - plane                                           *
*            CG(*)   - CG-object                                       *
*            RMN(3)  - min of CG-object scope                          *
*            RMX(3)  - max of CG-object scope                          *
*            NMAX    - max length of C array                           *
*                                                                      *
*     Output: C(*)   - resulting face                                  *
*             IREP   - reply (length of resulitng face C)              *
*                            (-1 if no space)                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      FACE(*),CG(*),RMN(3),RMX(3),C(*),D(8),ABCD(4)
*-
      IF (NMAX .LT. LCGFAC)     GOTO 998
      IREP   = 0
      C(KCGAF) = NMAX
      C(KCGNE) = 0.
      ABCD(1)= FACE(KCGAA)
      ABCD(2)= FACE(KCGBB)
      ABCD(3)= FACE(KCGCC)
      ABCD(4)= FACE(KCGDD)
*           M I N - M A X    T E S T
      D(1)   = ABCD(1)*RMN(1)+ABCD(2)*RMN(2)+ABCD(3)*RMN(3)+ABCD(4)
      D(2)   = ABCD(1)*RMX(1)+ABCD(2)*RMN(2)+ABCD(3)*RMN(3)+ABCD(4)
      D(3)   = ABCD(1)*RMN(1)+ABCD(2)*RMX(2)+ABCD(3)*RMN(3)+ABCD(4)
      D(4)   = ABCD(1)*RMX(1)+ABCD(2)*RMX(2)+ABCD(3)*RMN(3)+ABCD(4)
      D(5)   = ABCD(1)*RMN(1)+ABCD(2)*RMN(2)+ABCD(3)*RMX(3)+ABCD(4)
      D(6)   = ABCD(1)*RMX(1)+ABCD(2)*RMN(2)+ABCD(3)*RMX(3)+ABCD(4)
      D(7)   = ABCD(1)*RMN(1)+ABCD(2)*RMX(2)+ABCD(3)*RMX(3)+ABCD(4)
      D(8)   = ABCD(1)*RMX(1)+ABCD(2)*RMX(2)+ABCD(3)*RMX(3)+ABCD(4)
      NPOS   = 0
      NNEG   = 0
      DO 100 I=1,8
        IF (D(I) .GT. EEWOR)    NPOS = NPOS + 1
        IF (D(I) .LT.-EEWOR)    NNEG = NNEG + 1
  100   CONTINUE
      IF (NPOS.EQ.8 .OR. NNEG.EQ.8)     GOTO 999
*
      NFACE  = CG(KCGNF)
      JCG    = LCGHEA
      DO 200 NF=1,NFACE
        CALL CGBTFP(CG(JCG+1),ABCD,C)
        IF (C(KCGAF) .LT. 0.)   GOTO 998
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC + NEDGE*LCGEDG
  200   CONTINUE
*
**          S O R T   E G E S  &  D E L E T E   D O U B L E   E D G E S
*
      NEDGE  = C(KCGNE)
      CALL CGBSOR(NEDGE,C(LCGFAC+1))
      C(KCGNE) = NEDGE
      IF (NEDGE .EQ. 0)         GOTO 997
*      IF (NEDGE .LT. 3) PRINT *,' CGBSEC: NEDGE .LT. 3 - face ignored'
      IF (NEDGE .LT. 3)         GOTO 997
      C(KCGAA) = ABCD(1)
      C(KCGBB) = ABCD(2)
      C(KCGCC) = ABCD(3)
      C(KCGDD) = ABCD(4)
      J      = LCGFAC
      DO 300 NE=1,NEDGE
        C(J+KCGAE) = C(J+KCGAE) - 2.
        J    = J + LCGEDG
  300   CONTINUE
      IREP   = LCGFAC + LCGEDG*NEDGE
      GOTO 999
*
  997 IREP   = 0
      GOTO 999
  998 IREP   =-1
  999 RETURN
      END
