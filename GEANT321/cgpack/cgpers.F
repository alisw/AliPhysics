*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGPERS(CG)
************************************************************************
*                                                                      *
*     Name: CGPERS                                                     *
*     Author: S.Giani                            Date:    22.05.92     *
*                                                                      *
*     Function: Transform CG object according to  perspective view     *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input:  CG - CG-object                                           *
*     Output:                                                          *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
#include "geant321/gcdraw.inc"
      REAL      CG(*),GRAV(3)
      COMMON/PROSP/SVN(3)
*-
      IF (CG(KCGSIZ) .LE. 0.) GOTO 999
      LENG  = CG(KCGSIZ)
      CG(KCGSIZ) = 0.
      IF (LENG .LE. LCGHEA) GOTO 999
      NFACE  = CG(KCGNF)
      IF (NFACE .LE. 0) GOTO 999
      JCG    = LCGHEA
      DO 20  NF=1,NFACE
         JCGFAC = JCG
         NEDGE = CG(JCG+KCGNE)
         IF (NEDGE .LE. 0) GOTO 999
         A = 0.
         B = 0.
         C = 0.
         GRAV(1)= 0.
         GRAV(2)= 0.
         GRAV(3)= 0.
         JCG = JCG + LCGFAC
*             C O M P U T E   F A C E   A R E A
         DO 10  NE=1,NEDGE
*SG
*          PERSPECTIVE
*
            CALL UCTOH('PERS',IPERS,4,4)
            IF (IPRJ.EQ.IPERS) THEN
               XYZ3=CG(JCG+KCGX1)*SVN(1)+CG(JCG+KCGY1)*SVN(2)+ CG(JCG+
     +         KCGZ1)*SVN(3)
               IF (XYZ3.GE.DPERS) XYZ3=DPERS-0.1
               F=DPERS/(DPERS-XYZ3)
               CG(JCG+KCGX1)=CG(JCG+KCGX1)*F
               CG(JCG+KCGY1)=CG(JCG+KCGY1)*F
               CG(JCG+KCGZ1)=CG(JCG+KCGZ1)*F
               XYZ3=CG(JCG+KCGX2)*SVN(1)+CG(JCG+KCGY2)*SVN(2)+ CG(JCG+
     +         KCGZ2)*SVN(3)
               IF (XYZ3.GE.DPERS) XYZ3=DPERS-0.1
               F=DPERS/(DPERS-XYZ3)
               CG(JCG+KCGX2)=CG(JCG+KCGX2)*F
               CG(JCG+KCGY2)=CG(JCG+KCGY2)*F
               CG(JCG+KCGZ2)=CG(JCG+KCGZ2)*F
            ENDIF
*SG
            A = A + CG(JCG+KCGY1) * CG(JCG+KCGZ2) - CG(JCG+KCGY2) *
     +      CG(JCG+KCGZ1)
            B = B + CG(JCG+KCGZ1) * CG(JCG+KCGX2) - CG(JCG+KCGZ2) *
     +      CG(JCG+KCGX1)
            C = C + CG(JCG+KCGX1) * CG(JCG+KCGY2) - CG(JCG+KCGX2) *
     +      CG(JCG+KCGY1)
            GRAV(1)= GRAV(1) + CG(JCG+KCGX1) + CG(JCG+KCGX2)
            GRAV(2)= GRAV(2) + CG(JCG+KCGY1) + CG(JCG+KCGY2)
            GRAV(3)= GRAV(3) + CG(JCG+KCGZ1) + CG(JCG+KCGZ2)
            JCG = JCG + LCGEDG
   10    CONTINUE
         GRAV(1) = GRAV(1) / (2*NEDGE)
         GRAV(2) = GRAV(2) / (2*NEDGE)
         GRAV(3) = GRAV(3) / (2*NEDGE)
         AREA = SQRT(A*A + B*B + C*C)
*SG
*
*   Consider pyramid as a limit TRAP
*
         IF (AREA .LT. EEWOR)AREA=.1
*
         IF (AREA .LT. EEWOR) GOTO 999
*SG
         D =-(A*GRAV(1) + B*GRAV(2) + C*GRAV(3)) / AREA
         A = A / AREA
         B = B / AREA
         C = C / AREA
         CG(JCGFAC+KCGAA) = A
         CG(JCGFAC+KCGBB) = B
         CG(JCGFAC+KCGCC) = C
         CG(JCGFAC+KCGDD) = D
   20 CONTINUE
      CG(KCGSIZ) = LENG
  999 END
