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
      SUBROUTINE CGBOOL(IFA,IFB,A,B,NMAX,C)
************************************************************************
*                                                                      *
*     Name: CGBOOL                                                     *
*     Author: E. Chernyaev                       Date:    05.04.89     *
*                                                Revised:              *
*                                                                      *
*     Function: Make boolean operation with                            *
*                                                                      *
*     References: CGMNMX,CGCOPY,CGINVE,CGBSEA,CGBSEC,CGBINT,CGBSUB     *
*                                                                      *
*     Input: A(*) - 1-st CG-object                                     *
*            B(*) - 2-nd CG-object                                     *
*            NMAX - max length C array                                 *
*                                                                      *
*     Output: C(*) - resulting CG-object                               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL              A(*),B(*),C(*),AMN(3),AMX(3),BMN(3),BMX(3)
*-
      C(KCGSIZ)   = 0.
*           T E S T   P A R A M E T E R S   C O R R E C T N E S S
      CALL CGMNMX(A,AMN,AMX)
      CALL CGMNMX(B,BMN,BMX)
      IF (AMN(1) .GT. AMX(1))           GOTO 999
      IF (BMN(1) .GT. BMX(1))           GOTO 999
      IF (NMAX .LT. LCGHEA)             GOTO 998
      IF (A(KCGNF)  .EQ. 0.)            GOTO 110
      IF (B(KCGNF)  .EQ. 0.)            GOTO 120
      GOTO 150
*           O B J E C T   "A"   I S   E M P T Y
  110 IF (IFA.NE.0 .AND. IFB.EQ.0)      GOTO 999
      IF (IFA .EQ. 0)                   CALL CGCOPY(B,NMAX,C)
      IF (IFB .NE. 0)                   CALL CGCOPY(A,NMAX,C)
      IF (C(KCGSIZ) .LT. 0.)            GOTO 998
      GOTO 999
*           O B J E C T   "B"   I S   E M P T Y
  120 IF (IFB .EQ. 0)                   CALL CGCOPY(A,NMAX,C)
      IF (IFB .NE. 0)                   CALL CGCOPY(B,NMAX,C)
      IF (C(KCGSIZ) .LT. 0.)            GOTO 998
      GOTO 999
*
**          M I N - M A X   T E S T
*
  150 DO 200 I=1,3
        IF (AMN(I) .GT. BMX(I)+EEWOR)   GOTO 210
        IF (AMX(I) .LT. BMN(I)-EEWOR)   GOTO 210
  200   CONTINUE
      GOTO 300
*
**          N O    I N T E R S E C T I O N   O F   S C O P E S
*
  210 IF (IFA.NE.0 .AND. IFB.NE.0)      GOTO 220
****  12-MAR-1990
*
*      IF (IFA.NE.0 .AND. IFB.EQ.0)      GOTO 220
      IF (IFA.NE.0 .AND. IFB.EQ.0)      GOTO 230
*
      IF (IFA.EQ.0 .AND. IFB.EQ.0)      GOTO 240
*           N O   I N T E R S E C T I O N
  220 C(KCGSIZ) = LCGHEA
      C(KCGATT) = 0.
      C(KCGNF)  = 0.
      GOTO 999
*           N O   S U B T R A C T I O N
  230 CALL CGCOPY(A,NMAX,C)
      IF (C(KCGSIZ) .LT. 0.)            GOTO 998
      GOTO 999
*           S I M P L E   S U M   O F   O B J E C T S
  240 LTOTAL = A(KCGSIZ) + B(KCGSIZ) - LCGHEA
      IF (LTOTAL .GT. NMAX)             GOTO 998
      J      = A(KCGSIZ) - LCGHEA + 1
      CALL CGCOPY(B,NMAX,C(J))
      CALL CGCOPY(A,NMAX,C)
      C(KCGSIZ) = LTOTAL
      C(KCGNF)  = A(KCGNF) + B(KCGNF)
      GOTO 999
*
**          F I N D   I N T E R S E C T I O N   O F   F A C E S
**          W I T H   T H E   S A M E   P L A N E
*
  300 IF (IFA .NE. 0)           CALL CGINVE(A)
      IF (IFB .NE. 0)           CALL CGINVE(B)
      JA     = LCGHEA
      JC     = LCGHEA
      NC     = 0
      NFACE  = A(KCGNF)
      DO 330 NF=1,NFACE
        JB     = LCGHEA
        NB     = 0
  310   CALL CGBSEA(A(JA+1),B,NB,JB)
        IF (NB .EQ. 0)  GOTO 320
        CALL CGBINT(A(JA+1),B(JB+1),NMAX-JC,C(JC+1),IREP)
        IF (IREP .LT. 0)        GOTO 997
        IF (IREP .EQ. 0)        GOTO 310
        NC     = NC + 1
        JC     = JC + IREP
        GOTO 310
  320   NEDGE  = A(JA+KCGNE)
        JA     = JA + LCGFAC + NEDGE*LCGEDG
  330   CONTINUE
*
**          F I N D   P A R T   O F   1 - S T   O B J E C T   W H I C H
**          L I E   O U T S I D E   O F   2 - N D   O B J E C T
*
      J      = LCGHEA
      NFACE  = A(KCGNF)
      DO 400 NF=1,NFACE
        JCSAV  = JC
        CALL CGBSEC(A(J+1),B,BMN,BMX,NMAX-JCSAV,C(JCSAV+1),IREP)
        IF (IREP .LT. 0)        GOTO 997
        JC     = JCSAV + IREP
        CALL CGBSUB(IFB,A(J+1),C(JCSAV+1),NMAX-JC,C(JC+1),IREP)
        IF (IREP .LT. 0)        GOTO 997
        JC     = JCSAV + IREP
        IF (IREP .GT. 0)        NC = NC + 1
        NEDGE  = A(J+KCGNE)
        J      = J + LCGFAC + NEDGE*LCGEDG
  400   CONTINUE
*
**          F I N D   P A R T   O F   2 - N D   O B J E C T   W H I C H
**          L I E   O U T S I D E   O F   1 - S T   O B J E C T
*
      J      = LCGHEA
      NFACE  = B(KCGNF)
      DO 500 NF=1,NFACE
        JCSAV  = JC
        CALL CGBSEC(B(J+1),A,AMN,AMX,NMAX-JCSAV,C(JCSAV+1),IREP)
        IF (IREP .LT. 0)        GOTO 997
        JC     = JCSAV + IREP
        CALL CGBSUB(IFA,B(J+1),C(JCSAV+1),NMAX-JC,C(JC+1),IREP)
        IF (IREP .LT. 0)        GOTO 997
        JC     = JCSAV + IREP
        IF (IREP .GT. 0)        NC = NC + 1
        NEDGE  = B(J+KCGNE)
        J      = J + LCGFAC + NEDGE*LCGEDG
  500   CONTINUE
*
      C(KCGSIZ) = JC
      C(KCGATT) = 0.
      C(KCGNF)  = NC
      IF (IFA.NE.0 .OR. IFB.NE.0)       CALL CGINVE(C)
*
  997 IF (IFA .NE. 0)                   CALL CGINVE(A)
      IF (IFB .NE. 0)                   CALL CGINVE(B)
      IF (IREP .GE. 0)                  GOTO 999
  998 C(KCGSIZ) = -1.
  999 RETURN
      END
