*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:29  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.27  by  S.Giani
*-- Author :
      SUBROUTINE GDTR99(IVTREE)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Delete the tree structure build by GDTR0                 *
C.    *       on the view bank IVTREE                                  *
C.    *                                                                *
C.    *       IVTREE= view bank used to store tree structure (input)   *
C.    *                                                                *
C.    *    ==>Called by : GDFSPC, GDTREE                               *
C.    *       Author : P.Zanarini   *********                          *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gcdraw.inc"
C.
C.    ------------------------------------------------------------------
C.
      IF (NKVIEW.EQ.IVTREE) THEN
         Q(JDRAW+IVTREE)=1.
         CALL MZDROP(IXCONS,LQ(JDRAW-IVTREE),' ')
         NKVIEW=IVTREE-1
      ELSE
         Q(JDRAW+IVTREE)=1.
         CALL GDELET(IVTREE)
      ENDIF
      END
