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
      SUBROUTINE CGBMMT(A,B,IREP)
************************************************************************
*                                                                      *
*     Name: CGBMMT                                                     *
*     Author: E. Chernyaev                       Date:    23.10.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Mim-max test of two faces                              *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: A(*) - 1-st face                                          *
*            B(*) - 2-nd face                                          *
*                                                                      *
*     Output: IREP - reply (=0 if NO intersection)                     *
*                          (=1 if there is intersection)               *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
#include "geant321/cggpar.inc"
#include "geant321/cgdelt.inc"
      REAL      A(*),B(*),RMNA(3),RMXA(3),RMNB(3),RMXB(3)
      INTEGER   IND(6)
      DATA      IND/KCGX1,KCGY1,KCGZ1,KCGX2,KCGY2,KCGZ2/
*-
      JJ     =  LCGFAC
      DO 100 I=1,3
        RMNA(I)= A(JJ+IND(I))
        RMXA(I)= A(JJ+IND(I))
        RMNB(I)= B(JJ+IND(I))
        RMXB(I)= B(JJ+IND(I))
  100   CONTINUE
*           C O M P U T E   M I N - M A X   O F   1 - S T   F A C E
      NEDGE  = A(KCGNE)
      JJ     = LCGFAC
      DO 220 NE=1,NEDGE
        DO 210 J=1,3
          I      = J
          IF (J .GT. 3) I = J - 3
          IF (A(JJ+IND(J)) .LT. RMNA(I)) RMNA(I) = A(JJ+IND(J))
          IF (A(JJ+IND(J)) .GT. RMXA(I)) RMXA(I) = A(JJ+IND(J))
  210     CONTINUE
        JJ     = JJ + LCGEDG
  220   CONTINUE
*           C O M P U T E   M I N - M A X   O F   2 - N D   F A C E
      NEDGE  = B(KCGNE)
      JJ     = LCGFAC
      DO 320 NE=1,NEDGE
        DO 310 J=1,3
          I      = J
          IF (J .GT. 3) I = J - 3
          IF (B(JJ+IND(J)) .LT. RMNB(I)) RMNB(I) = B(JJ+IND(J))
          IF (B(JJ+IND(J)) .GT. RMXB(I)) RMXB(I) = B(JJ+IND(J))
  310     CONTINUE
        JJ     = JJ + LCGEDG
  320   CONTINUE
*           M I N - M A X   T E S T
      IREP   = 0
      IF (RMNA(1) .GT. RMXB(1)+EEWOR)   GOTO 999
      IF (RMXA(1) .LT. RMNB(1)-EEWOR)   GOTO 999
      IF (RMNA(2) .GT. RMXB(2)+EEWOR)   GOTO 999
      IF (RMXA(2) .LT. RMNB(2)-EEWOR)   GOTO 999
      IF (RMNA(3) .GT. RMXB(3)+EEWOR)   GOTO 999
      IF (RMXA(3) .LT. RMNB(3)-EEWOR)   GOTO 999
      IREP   = 1
*
  999 RETURN
      END
