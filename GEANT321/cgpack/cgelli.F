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
      SUBROUTINE CGELLI(RX,RY,RZ,KA,KB,LCG,CG)
************************************************************************
*                                                                      *
*     Name: CGELLI                                                     *
*     Author: E. Chernyaev                       Date:    24.01.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Create CG-object for ELLIPSOID                         *
*                                                                      *
*     References: CGSIZE, CGSNOR                                       *
*                                                                      *
*     Input:  RX - 1-st radius                                         *
*             RY - 2-nd radius                                         *
*             RZ - 3-rd radius                                         *
*             KA - number of latitude step                             *
*             KB - number of longitude step                            *
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
      REAL      CG(*)
      INTEGER   NFAC(2),NEDG(2)
*-
      CG(KCGSIZ)  = 0.
*           T E S T   P A R A M E T E R S   C O R R E C T N E S S
      IF (RX .LE. 0.)   GOTO 999
      IF (RY .LE. 0.)   GOTO 999
      IF (RZ .LE. 0.)   GOTO 999
      IF (KA .LT. 3)    GOTO 999
      IF (KB .LT. 2)    GOTO 999
*           C O M P U T E   S I Z E   O F   C G - O B J E C T
      NFATYP = 2
      NFAC(1)= 2 * KA
      NEDG(1)= 3
      NFAC(2)= (KB - 2) * KA
      NEDG(2)= 4
      CG(KCGSIZ)  = CGSIZE(LCG,NFATYP,NFAC,NEDG)
      IF (CG(KCGSIZ) .LE. 0.)        GOTO 999
*
**          C R E A T E   C G - O B J E C T   F O R   E L L I P S O I D
*
      CG(KCGATT) = 0.
      CG(KCGNF)  = NFAC(1) + NFAC(2)
      PI     = ATAN(1.) * 4.
      DA     = (PI+PI) / KA
      DB     = PI / KB
      JCOSB  = CG(KCGSIZ) - 2*(KB+KA) - 4
      JSINB  = JCOSB + KB + 1
      JCOSA  = JSINB + KB + 1
      JSINA  = JCOSA + KA + 1
*           P R E P A R E   W O R K   T A B L E
      A      = 0.
      B      = 0.
      JA     = 0
      JB     = 0
  100 CG(JCOSA+JA)     = COS(A)
      CG(JCOSA+KA-JA)  = CG(JCOSA+JA)
      CG(JSINA+JA)     = SIN(A)
      CG(JSINA+KA-JA)  =-CG(JSINA+JA)
      JA     = JA + 1
      A      = A + DA
      IF (KA-JA-JA) 200,110,100
  110 CG(JCOSA+JA)     =-1.
      CG(JSINA+JA)     = 0.
*
  200 CG(JCOSB+JB)    = COS(B)
      CG(JCOSB+KB-JB) =-CG(JCOSB+JB)
      CG(JSINB+JB)    = SIN(B)
      CG(JSINB+KB-JB) = CG(JSINB+JB)
      JB     = JB + 1
      B      = B + DB
      IF (KB-JB-JB) 300,210,200
  210 CG(JCOSB+JB)    = 0.
      CG(JSINB+JB)    = 1.
*           C R E A T E   U P P E R   H A L F   O F   E L L I P S O I D
  300 JCG    = LCGHEA
      NB     = KB - KB/2
      DO 500 JB=1,NB
        Z1     = RZ*CG(JCOSB+JB-1)
        Z2     = RZ*CG(JCOSB+JB)
        X4     = RX*CG(JSINB+JB-1)
        X3     = RX*CG(JSINB+JB)
        Y4     = 0.
        Y3     = 0.
        DO 400 JA=1,KA
          CG(JCG+KCGAF) = 0.
          JCGNE  = JCG + KCGNE
          JCG    = JCG + LCGFAC
          X1     = X4
          X2     = X3
          X3     = RX*CG(JSINB+JB)   * CG(JCOSA+JA)
          X4     = RX*CG(JSINB+JB-1) * CG(JCOSA+JA)
          Y1     = Y4
          Y2     = Y3
          Y3     = RY*CG(JSINB+JB)   * CG(JSINA+JA)
          Y4     = RY*CG(JSINB+JB-1) * CG(JSINA+JA)
          CG(JCG+KCGAE) =-1.
          CG(JCG+KCGX1) = X1
          CG(JCG+KCGY1) = Y1
          CG(JCG+KCGZ1) = Z1
          CG(JCG+KCGX2) = X2
          CG(JCG+KCGY2) = Y2
          CG(JCG+KCGZ2) = Z2
*
          CG(JCG+LCGEDG+KCGAE) =-1.
          CG(JCG+LCGEDG+KCGX1) = X2
          CG(JCG+LCGEDG+KCGY1) = Y2
          CG(JCG+LCGEDG+KCGZ1) = Z2
          CG(JCG+LCGEDG+KCGX2) = X3
          CG(JCG+LCGEDG+KCGY2) = Y3
          CG(JCG+LCGEDG+KCGZ2) = Z2
*
          CG(JCG+LCGEDG+LCGEDG+KCGAE) =-1.
          CG(JCG+LCGEDG+LCGEDG+KCGX1) = X3
          CG(JCG+LCGEDG+LCGEDG+KCGY1) = Y3
          CG(JCG+LCGEDG+LCGEDG+KCGZ1) = Z2
          CG(JCG+LCGEDG+LCGEDG+KCGX2) = X4
          CG(JCG+LCGEDG+LCGEDG+KCGY2) = Y4
          CG(JCG+LCGEDG+LCGEDG+KCGZ2) = Z1
          NE     = 3
          IF(X1.EQ.X4 .AND. Y1.EQ.Y4) GOTO 350
          NE     = 4
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGAE) =-1.
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGX1) = X4
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGY1) = Y4
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGZ1) = Z1
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGX2) = X1
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGY2) = Y1
          CG(JCG+LCGEDG+LCGEDG+LCGEDG+KCGZ2) = Z1
  350     JCG    = JCG + NE*LCGEDG
          CG(JCGNE) = NE
  400     CONTINUE
  500   CONTINUE
*           C R E A T E   L O W E R   H A L F   O F   E L L I P S O I D
      JSTOP  = JCG
      JCGUP  = LCGHEA
      JCGLOW = CG(KCGSIZ)
  600 NE     = CG(JCGUP+KCGNE)
      JCGLOW = JCGLOW - LCGFAC - NE*LCGEDG
      IF (JCGLOW .LT. JSTOP)    GOTO 999
      CG(JCGLOW+KCGAF) = CG(JCGUP+KCGAF)
      CG(JCGLOW+KCGNE) = NE
      JL     = JCGLOW + LCGFAC
      JU     = JCGUP  + LCGFAC
*
      CG(JL+KCGAE) = CG(JU+KCGAE)
      CG(JL+KCGX1) = CG(JU+KCGX2)
      CG(JL+KCGY1) = CG(JU+KCGY2)
      CG(JL+KCGZ1) =-CG(JU+KCGZ2)
      CG(JL+KCGX2) = CG(JU+KCGX1)
      CG(JL+KCGY2) = CG(JU+KCGY1)
      CG(JL+KCGZ2) =-CG(JU+KCGZ1)
      JU = JU + (NE-1)*LCGEDG
*
      CG(JL+LCGEDG+KCGAE) = CG(JU+KCGAE)
      CG(JL+LCGEDG+KCGX1) = CG(JU+KCGX2)
      CG(JL+LCGEDG+KCGY1) = CG(JU+KCGY2)
      CG(JL+LCGEDG+KCGZ1) =-CG(JU+KCGZ2)
      CG(JL+LCGEDG+KCGX2) = CG(JU+KCGX1)
      CG(JL+LCGEDG+KCGY2) = CG(JU+KCGY1)
      CG(JL+LCGEDG+KCGZ2) =-CG(JU+KCGZ1)
*
      CG(JL+LCGEDG+LCGEDG+KCGAE) = CG(JU-LCGEDG+KCGAE)
      CG(JL+LCGEDG+LCGEDG+KCGX1) = CG(JU-LCGEDG+KCGX2)
      CG(JL+LCGEDG+LCGEDG+KCGY1) = CG(JU-LCGEDG+KCGY2)
      CG(JL+LCGEDG+LCGEDG+KCGZ1) =-CG(JU-LCGEDG+KCGZ2)
      CG(JL+LCGEDG+LCGEDG+KCGX2) = CG(JU-LCGEDG+KCGX1)
      CG(JL+LCGEDG+LCGEDG+KCGY2) = CG(JU-LCGEDG+KCGY1)
      CG(JL+LCGEDG+LCGEDG+KCGZ2) =-CG(JU-LCGEDG+KCGZ1)
      IF (NE .EQ. 3)    GOTO 700
*
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGAE) = CG(JU-LCGEDG-LCGEDG+KCGAE)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGX1) = CG(JU-LCGEDG-LCGEDG+KCGX2)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGY1) = CG(JU-LCGEDG-LCGEDG+KCGY2)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGZ1) =-CG(JU-LCGEDG-LCGEDG+KCGZ2)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGX2) = CG(JU-LCGEDG-LCGEDG+KCGX1)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGY2) = CG(JU-LCGEDG-LCGEDG+KCGY1)
      CG(JL+LCGEDG+LCGEDG+LCGEDG+KCGZ2) =-CG(JU-LCGEDG-LCGEDG+KCGZ1)
  700 JCGUP  = JCGUP + LCGFAC + NE*LCGEDG
      GOTO 600
*
  999 CALL CGSNOR(CG)
      RETURN
      END
