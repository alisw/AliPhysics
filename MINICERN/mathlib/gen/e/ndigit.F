*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:26  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
      FUNCTION NDIGIT(U,V)

#include "gen/imp64.inc"

************************************************************************
*   NORBAS, VERSION: 15.03.1993
************************************************************************
*
*   NDIGIT  COMPUTES THE NUMBER OF DECIMAL DIGITS AT THE BEGINNING OF  U
*   WHICH ARE AGREE WITH CORRESPONDING DIGITS OF  V
*
************************************************************************

      PARAMETER (SMALL = 1D-30 , TWO = 2)
      IF(ABS(U) .LE. SMALL .AND. ABS(V) .LE. SMALL) THEN
       NDIGIT=16
      ELSE IF(U*V .LT. 0) THEN
       NDIGIT=0
      ELSE
       NDIGIT=MIN(16,ABS(INT(LOG10(MAX(SMALL,ABS((U+V))))) -
     1                   INT(LOG10(MAX(SMALL,TWO*ABS(U-V))))))
      ENDIF
      RETURN
      END

