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
      SUBROUTINE CGBOX (XYZ,N1,N2,LCG,CG)
************************************************************************
*                                                                      *
*     Name: CGBOX                                                      *
*     Author: E. Chernyaev                       Date:    12.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Create CG-object for box                               *
*                                                                      *
*     References: CGSIZE, CGSNOR                                       *
*                                                                      *
*     Input: XYZ(3,N1+N2) - nodes                                      *
*                N1 - number of nodes in main base                     *
*                N2 - number of nodes in 2-nd base                     *
*               LCG - max-size of CG-object                            *
*                                                                      *
*     Output: CG - CG-object                                           *
*                  CG(1) - length of CG-object                         *
*                        = 0 if error in parameters                    *
*                        < 0 if no space                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
      REAL      CG(*),XYZ(3,*)
      INTEGER   NFAC(2),NEDG(2)
*-
      CG(KCGSIZ)  = 0.
*           T E S T   P A R A M E T E R S   C O R R E C T N E S S
      IF (N1 .LT. 3)                    GOTO 999
      IF (N2.NE.N1 .AND. N2.NE.1)       GOTO 999
*           C O M P U T E   S I Z E   O F   C G - O B J E C T
      NFATYP = 2
      NEDG(1)= N1
      IF (N2 .NE. 1)                    NEDG(2) = 4
      IF (N2 .EQ. 1)                    NEDG(2) = 3
      IF (N2 .NE. 1)                    NFAC(1) = 2
      IF (N2 .EQ. 1)                    NFAC(1) = 1
      NFAC(2)= N1
      CG(KCGSIZ)  = CGSIZE(LCG,NFATYP,NFAC,NEDG)
      IF (CG(KCGSIZ) .LE. 0.)        GOTO 999
*
**          C R E A T E   C G - O B J E C T
*
      CG(KCGATT) = 0.
      CG(KCGNF)  = NFAC(1) + NFAC(2)
      JCG        = LCGHEA
*           C R E A T E   M A I N   B A S E
      CG(JCG + KCGAF) = 0.
      CG(JCG + KCGNE) = N1
      JCG    = JCG + LCGFAC
      DO 100 I=1,N1
        I1     = I
        I2     = I + 1
        IF (I .EQ. N1)          I2 = 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,I1)
        CG(JCG+KCGY1) = XYZ(2,I1)
        CG(JCG+KCGZ1) = XYZ(3,I1)
        CG(JCG+KCGX2) = XYZ(1,I2)
        CG(JCG+KCGY2) = XYZ(2,I2)
        CG(JCG+KCGZ2) = XYZ(3,I2)
        JCG    = JCG + LCGEDG
  100   CONTINUE
*           C R E A T E   2 - N D   B A S E
      IF (N2 .EQ. 1)            GOTO 300
      CG(JCG + KCGAF) = 0.
      CG(JCG + KCGNE) = N2
      JCG    = JCG + LCGFAC
      DO 200 I=1,N2
        I1     = I
        I2     = I + 1
        IF (I .EQ. N2)          I2 = 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,I2+N1)
        CG(JCG+KCGY1) = XYZ(2,I2+N1)
        CG(JCG+KCGZ1) = XYZ(3,I2+N1)
        CG(JCG+KCGX2) = XYZ(1,I1+N1)
        CG(JCG+KCGY2) = XYZ(2,I1+N1)
        CG(JCG+KCGZ2) = XYZ(3,I1+N1)
        JCG    = JCG + LCGEDG
  200   CONTINUE
*           C R E A T E   S I D E   F A C E S
  300 DO 350 I=1,N1
        I1     = I
        I2     = I + 1
        IF (I .EQ. N1)          I2 = 1
        CG(JCG + KCGAF) = 0.
        CG(JCG + KCGNE) = NEDG(2)
        JCG    = JCG + LCGFAC
*
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,I2)
        CG(JCG+KCGY1) = XYZ(2,I2)
        CG(JCG+KCGZ1) = XYZ(3,I2)
        CG(JCG+KCGX2) = XYZ(1,I1)
        CG(JCG+KCGY2) = XYZ(2,I1)
        CG(JCG+KCGZ2) = XYZ(3,I1)
        JCG    = JCG + LCGEDG
*
        IF (N2 .NE. 1)          K = N1 + I1
        IF (N2 .EQ. 1)          K = N1 + 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,I1)
        CG(JCG+KCGY1) = XYZ(2,I1)
        CG(JCG+KCGZ1) = XYZ(3,I1)
        CG(JCG+KCGX2) = XYZ(1,K)
        CG(JCG+KCGY2) = XYZ(2,K)
        CG(JCG+KCGZ2) = XYZ(3,K)
        JCG    = JCG + LCGEDG
*
        IF (N2 .EQ. 1)          GOTO 310
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,I1+N1)
        CG(JCG+KCGY1) = XYZ(2,I1+N1)
        CG(JCG+KCGZ1) = XYZ(3,I1+N1)
        CG(JCG+KCGX2) = XYZ(1,I2+N1)
        CG(JCG+KCGY2) = XYZ(2,I2+N1)
        CG(JCG+KCGZ2) = XYZ(3,I2+N1)
        JCG    = JCG + LCGEDG
*
  310   IF (N2 .NE. 1)          K = N1 + I2
        IF (N2 .EQ. 1)          K = N1 + 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = XYZ(1,K)
        CG(JCG+KCGY1) = XYZ(2,K)
        CG(JCG+KCGZ1) = XYZ(3,K)
        CG(JCG+KCGX2) = XYZ(1,I2)
        CG(JCG+KCGY2) = XYZ(2,I2)
        CG(JCG+KCGZ2) = XYZ(3,I2)
        JCG    = JCG + LCGEDG
  350   CONTINUE
*           S E T   N O R M A L E S
      CALL CGSNOR(CG)
      IF (CG(1) .EQ. 0.) PRINT *, ' Problem in CGBOX'
      JCG    = LCGHEA
      A      = CG(JCG + KCGAA)
      B      = CG(JCG + KCGBB)
      C      = CG(JCG + KCGCC)
      D      = CG(JCG + KCGDD)
      DIST   = A*XYZ(1,N1+1) + B*XYZ(2,N1+1) + C*XYZ(3,N1+1) + D
      IF (DIST .LT. 0.)         GOTO 999
*           S E T   I N V E R S E   O R D E R   O F   E D G E S
      NFACE  = CG(KCGNF)
      DO 500 NF=1,NFACE
        CG(JCG + KCGAA) =-CG(JCG + KCGAA)
        CG(JCG + KCGBB) =-CG(JCG + KCGBB)
        CG(JCG + KCGCC) =-CG(JCG + KCGCC)
        CG(JCG + KCGDD) =-CG(JCG + KCGDD)
        NEDGE  = CG(JCG+KCGNE)
        JCG    = JCG + LCGFAC
        DO 400 NE=1,NEDGE
          X      = CG(JCG+KCGX1)
          Y      = CG(JCG+KCGY1)
          Z      = CG(JCG+KCGZ1)
          CG(JCG + KCGX1) = CG(JCG + KCGX2)
          CG(JCG + KCGY1) = CG(JCG + KCGY2)
          CG(JCG + KCGZ1) = CG(JCG + KCGZ2)
          CG(JCG + KCGX2) = X
          CG(JCG + KCGY2) = Y
          CG(JCG + KCGZ2) = Z
          JCG    = JCG + LCGEDG
  400     CONTINUE
  500   CONTINUE
*
  999 RETURN
      END
