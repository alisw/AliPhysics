*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:57  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.43  by  S.Giani
*-- Author :
*$ CREATE HEVHIN.FOR
*COPY HEVHIN
*
*=== hevhin ===========================================================*
*
      SUBROUTINE HEVHIN ( NHAD, KPROJ, KTARG, PPROJ, EPROJ, UMO )
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*
*----------------------------------------------------------------------*
*                                                                      *
*     Created on 01 october 1991   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 01-oct-91     by    Alfredo Ferrari               *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*
#include "geant321/finlsp.inc"
#include "geant321/hadpar.inc"
#include "geant321/part.inc"
      PARAMETER ( CXDUM = 0.D+00 )
      PARAMETER ( CYDUM = 0.D+00 )
      PARAMETER ( CZDUM = 1.D+00 )
*
      ELAB = EPROJ
      PLAB = PPROJ
      N    = KPROJ
      ITTA = KTARG
      CALL HADRIV ( N, PLAB, ELAB, CXDUM, CYDUM, CZDUM, ITTA )
      NHAD = IR
      DO 1000 I = 1, IR
         HEPH (I) = ELR (I)
         PXH  (I) = PLR (I) * CXR (I)
         PYH  (I) = PLR (I) * CYR (I)
         PZH  (I) = PLR (I) * CZR (I)
         NREH (I) = ITR (I)
         ICHH (I) = ICH  (NREH(I))
         IBARH(I) = IBAR (NREH(I))
         AMH  (I) = AM   (NREH(I))
         ANH  (I) = ANAME(NREH(I))
 1000 CONTINUE
*=== End of subroutine Hevhin =========================================*
      RETURN
      END
