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
      SUBROUTINE CGBRIK(DX,DY,DZ,LCG,CG)
************************************************************************
*                                                                      *
*     Name: CGBRIK                                                     *
*     Author: E. Chernyaev                       Date:    24.03.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Create CG-object for BRICK                             *
*                                                                      *
*     References: CGSIZE                                               *
*                                                                      *
*     Input:  DX - X-size                                              *
*             DY - Y-size                                              *
*             DZ - Z-size                                              *
*            LCG - max-size of CG-object                               *
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
      REAL      CG(*),XYZ(3,8)
      INTEGER   IORDER(4,6)
      SAVE XYZ,IORDER
      DATA      XYZ/24*0./
      DATA      IORDER/1,2,3,4,
     +                 1,5,6,2, 2,6,7,3, 3,7,8,4, 4,8,5,1,
     +                 8,7,6,5/
*-
      CG(KCGSIZ)  = 0.
*           T E S T   P A R A M E T E R S   C O R R E C T N E S S
      IF (DX.LE.0. .OR. DY.LE.0. .OR. DZ.LE.0.) GOTO 999
*           C O M P U T E   S I Z E   O F   C G - O B J E C T
      NFATYP   = 1
      NFAC     = 6
      NEDG     = 4
      CG(KCGSIZ)  = CGSIZE(LCG,NFATYP,NFAC,NEDG)
      IF (CG(KCGSIZ) .LE. 0.)        GOTO 999
*           P R E P A R E   N O D E S   C O O R D I N A T E S
      XYZ(2,2) =  DY
      XYZ(1,3) =  DX
      XYZ(2,3) =  DY
      XYZ(1,4) =  DX
      DO 100 I=1,4
        XYZ(1,I+4) = XYZ(1,I)
        XYZ(2,I+4) = XYZ(2,I)
        XYZ(3,I+4) = DZ
  100   CONTINUE
*           C R E A T E   C G - O B J E C T   F O R   B R I K
      CG(KCGATT) = 0.
      CG(KCGNF)  = 6.
      JCG    = LCGHEA
      DO 300 NF=1,6
        CG(JCG+KCGAF) = 0.
        X1     = XYZ(1,IORDER(2,NF))-XYZ(1,IORDER(1,NF))
        Y1     = XYZ(2,IORDER(2,NF))-XYZ(2,IORDER(1,NF))
        Z1     = XYZ(3,IORDER(2,NF))-XYZ(3,IORDER(1,NF))
        X2     = XYZ(1,IORDER(3,NF))-XYZ(1,IORDER(2,NF))
        Y2     = XYZ(2,IORDER(3,NF))-XYZ(2,IORDER(2,NF))
        Z2     = XYZ(3,IORDER(3,NF))-XYZ(3,IORDER(2,NF))
        A      = Y1*Z2 - Y2*Z1
        IF (A .GT. 0.)  A = 1.
        IF (A .LT. 0.)  A =-1.
        B      = X2*Z1 - X1*Z2
        IF (B .GT. 0.)  B = 1.
        IF (B .LT. 0.)  B =-1.
        C      = X1*Y2 - X2*Y1
        IF (C .GT. 0.)  C = 1.
        IF (C .LT. 0.)  C =-1.
        D      =-(A*XYZ(1,IORDER(1,NF)) +
     +            B*XYZ(2,IORDER(1,NF)) +
     +            C*XYZ(3,IORDER(1,NF)) )
        CG(JCG+KCGAA) = A
        CG(JCG+KCGBB) = B
        CG(JCG+KCGCC) = C
        CG(JCG+KCGDD) = D
        CG(JCG+KCGNE) = 4.
        JCG    = JCG + LCGFAC
        DO 200 NE=1,4
          NE2    = NE + 1
          IF (NE .EQ. 4) NE2 = 1
          CG(JCG+KCGAE) = 0.
          CG(JCG+KCGX1) = XYZ(1,IORDER(NE,NF))
          CG(JCG+KCGY1) = XYZ(2,IORDER(NE,NF))
          CG(JCG+KCGZ1) = XYZ(3,IORDER(NE,NF))
          CG(JCG+KCGX2) = XYZ(1,IORDER(NE2,NF))
          CG(JCG+KCGY2) = XYZ(2,IORDER(NE2,NF))
          CG(JCG+KCGZ2) = XYZ(3,IORDER(NE2,NF))
          JCG    = JCG + LCGEDG
  200     CONTINUE
  300   CONTINUE
*
  999 RETURN
      END
