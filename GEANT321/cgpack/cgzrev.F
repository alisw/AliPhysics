*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
      SUBROUTINE CGZREV(RZ,A1,A2,NA,LCG,CG)
************************************************************************
*                                                                      *
*     Name: CGZREV                                                     *
*     Author: E. Chernyaev                       Date:    05.02.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Create CG-object by revolution around Z-axis           *
*                                                                      *
*     References: CGSIZE, CGSAAN, CGZRE                                *
*                                                                      *
*     Input: RZ(2,4) - 4 node contour (1-st must be left lower node)   *
*                A1 - initial angle                                    *
*                A2 - end angle                                        *
*                NA - number of steps on angle                         *
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
#include "geant321/cgdelt.inc"
#include "geant321/cgcaan.inc"
      REAL      RZ(2,4),CG(*)
#if !defined(CERNLIB_SINGLE)
      DOUBLE PRECISION  SINE,COSE
#endif
      INTEGER   NFAC(4),NEDG(4)
*-
      CG(KCGSIZ)  = 0.
*           T E S T   P A R A M E T E R S   C O R R E C T N E S S
      DO 100 I=1,4
        IF (RZ(1,I) .LT. 0.)    GOTO 999
        K      = I + 1
        IF (I .EQ. 4)           K = 1
        RLENG  = ABS(RZ(1,I)-RZ(1,K))
        ZLENG  = ABS(RZ(2,I)-RZ(2,K))
        IF (RLENG.LT.EEWOR .AND. ZLENG.LT.EEWOR)        GOTO 999
  100   CONTINUE
      CALL CGSAAN(A1,A2,NA,NA,IREP)
*           C O M P U T E   S I Z E   O F   C G - O B J E C T
      NFATYP = 4
      NEDG(1)= 3
      NEDG(2)= 4
      IF (IFULL .NE. 0)         NEDG(3)= NA
      IF (IFULL .EQ. 0)         NEDG(3)= NA + 2
      IF (IFULL .NE. 0)         NEDG(4)= NA + NA
      IF (IFULL .EQ. 0)         NEDG(4)= NA + NA + 2
      NFAC(1)= 0
      IF (IFULL .NE. 0)         NFAC(2)= 0
      IF (IFULL .EQ. 0)         NFAC(2)= 2
      NFAC(3)= 0
      NFAC(4)= 0
*
      DO 150 I=1,4
        K      = I + 1
        IF (I .EQ. 4)           K = 1
        IF (RZ(1,I).LT.EEWOR .AND. RZ(1,K).LT.EEWOR)    GOTO 150
        IF (ABS(RZ(2,I)-RZ(2,K)) .LT. EEWOR)            GOTO 110
        J      = 2
        IF (RZ(1,I) .LT. EEWOR)    J = J - 1
        IF (RZ(1,K) .LT. EEWOR)    J = J - 1
        NFAC(J) = NFAC(J) + NA
        GOTO 150
  110   IF (ABS(RZ(1,I)-RZ(1,K)) .LT. EEWOR)            GOTO 150
        J      = 4
        IF (RZ(1,I) .LT. EEWOR)    J = J - 1
        IF (RZ(1,K) .LT. EEWOR)    J = J - 1
        NFAC(J) = NFAC(J) + 1
  150   CONTINUE
      CG(KCGSIZ)  = CGSIZE(LCG,NFATYP,NFAC,NEDG)
      IF (CG(KCGSIZ) .LE. 0.)   GOTO 999
*
**          C R E A T E   C G - O B J E C T
*
      CG(KCGATT) = 0.
      CG(KCGNF)  = NFAC(1) + NFAC(2) + NFAC(3) + NFAC(4)
      JCG        = LCGHEA
      ATREDG     =-1.
      XYHA(2,1)  = 0.
      XYHA(2,2)  = 0.
      XYHA(4,1)  = 0.
      XYHA(4,2)  = 0.
      DO 200 I=1,4
        K      = I + 1
        IF (I .EQ. 4)           K = 1
        XYHA(1,1) = RZ(1,I)
        XYHA(3,1) = RZ(2,I)
        XYHA(1,2) = RZ(1,K)
        XYHA(3,2) = RZ(2,K)
        CALL CGZRE(2,CG(JCG+1),J)
        JCG    = JCG + J
  200   CONTINUE
      IF (IFULL .NE. 0.)        GOTO 999
*           C R E A T E   S I D E   F A C E S
      CG(JCG+KCGAF) = 0.
      CG(JCG+KCGAA) = SINI
      CG(JCG+KCGBB) =-COSI
      CG(JCG+KCGCC) = 0.
      CG(JCG+KCGDD) = 0.
      CG(JCG+KCGNE) = 4.
      JCG    = JCG + LCGFAC
      DO 300 I=1,4
        K      = I + 1
        IF (I .EQ. 4)           K = 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = RZ(1,I)*COSI
        CG(JCG+KCGY1) = RZ(1,I)*SINI
        CG(JCG+KCGZ1) = RZ(2,I)
        CG(JCG+KCGX2) = RZ(1,K)*COSI
        CG(JCG+KCGY2) = RZ(1,K)*SINI
        CG(JCG+KCGZ2) = RZ(2,K)
        JCG    = JCG + LCGEDG
  300   CONTINUE
*
      CG(JCG+KCGAF) = 0.
      COSE   = COS(A2*RAD)
      SINE   = SIN(A2*RAD)
      CG(JCG+KCGAA) =-SINE
      CG(JCG+KCGBB) = COSE
      CG(JCG+KCGCC) = 0.
      CG(JCG+KCGDD) = 0.
      CG(JCG+KCGNE) = 4.
      JCG    = JCG + LCGFAC
      DO 400 I=1,4
        K      = I + 1
        IF (I .EQ. 4)           K = 1
        CG(JCG+KCGAE) = 0.
        CG(JCG+KCGX1) = RZ(1,K)*COSE
        CG(JCG+KCGY1) = RZ(1,K)*SINE
        CG(JCG+KCGZ1) = RZ(2,K)
        CG(JCG+KCGX2) = RZ(1,I)*COSE
        CG(JCG+KCGY2) = RZ(1,I)*SINE
        CG(JCG+KCGZ2) = RZ(2,I)
        JCG    = JCG + LCGEDG
  400   CONTINUE
*
  999 RETURN
      END
