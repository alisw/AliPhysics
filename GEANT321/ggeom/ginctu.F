*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:50  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.29  by  S.Giani
*-- Author :
      SUBROUTINE GINCTU( X     , PAR   , IYES  )
*
***********************************************************************
*.
*...GINCTU tests whether a point inside a cut tube or not.
*.
*. SUBROUTINE GINCTU( X     , PAR   , IYES  )
*. -----------------
*.
*. Routine GINCTU tests whether a point is inside a cur tue or not.
*.
*. INPUT     : X      vector    real length at least 3, x,y,z of
*.           :                  point in local cut tube coords.
*.           : PAR    vector    real length 11 parameters of cut
*.           :                  tube.
*. OUTPUT    : IYES   integer   1 if point inside, 0 if not.
*.
*. SEQUENCES : GCONSP
*. CALLED    : GINME
*.
*.
*. REPORT CONDITIONS
*.
*. Non Stand.: none.
*.
*. AUTHOR    : A.C.McPherson.
*. VERSION   : 1.00
*. CREATED   :  6-May-88
*. LAST MOD  :  6-May-88
*.
*. MODIFICATION LOG.
*.
*.*********************************************************************
*.
#include "geant321/gconsp.inc"
      REAL X(*), PAR(*)
*.
*          --------------------
*
      IYES = 0
      ZINL = X(1)*PAR(6) +X(2)*PAR(7) +( X(3)+PAR(3) )*PAR(8)
      IF( ZINL .LT. 0.0 ) THEN
        ZINH = X(1)*PAR(9) +X(2)*PAR(10) +( X(3)-PAR(3) )*PAR(11)
        IF( ZINH .LT. 0.0 ) THEN
          R2 = X(1)*X(1)+X(2)*X(2)
          IF( R2 .GT. PAR(1)*PAR(1) .AND. R2 .LT. PAR(2)*PAR(2) ) THEN
            IF( PAR(4) .EQ. 0 .AND. PAR(5) .EQ. 360.0 ) THEN
              IYES = 1
            ELSEIF( R2 .GT. 0.0 ) THEN
              PHI = ATAN2( X(2), X(1) )/DEGRAD
              DPHI = PHI - PAR(4)
              IF( DPHI .LT. 0.0 ) DPHI = DPHI+360.0
              IF( DPHI .GT. 360.0 ) DPHI = DPHI-360.0
              DPT = PAR(5) - PAR(4)
              IF( DPT .LT. 0.0 ) DPT = DPT+360.0
              IF( DPT .GT. 360.0 ) DPT = DPT-360.0
              IF( DPHI .LE. DPT ) IYES = 1
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      END
