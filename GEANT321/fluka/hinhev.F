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
*$ CREATE HINHEV.FOR
*COPY HINHEV
*
*=== hinhev ===========================================================*
*
      SUBROUTINE HINHEV ( N, PLAB, ELAB, ITTA )
 
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
*
      EPROJ = ELAB
      PPROJ = PLAB
      KPROJ = N
      KTARG = ITTA
      AMTAR = AM (KTARG)
      UMO   = SQRT ( ( EPROJ + AMTAR - PPROJ ) * ( EPROJ + AMTAR
     &      + PPROJ ) )
      NHAD  = 0
      CALL HADEVV ( NHAD, KPROJ, KTARG, PPROJ, EPROJ, UMO )
      IR = NHAD
      DO 1000 I = 1, NHAD
         ELR (I) = HEPH (I)
         PLR (I) = SQRT ( ( HEPH (I) - AMH (I) ) * ( HEPH (I) + AMH (I)
     &           ) )
         CXR (I) = PXH (I) / PLR (I)
         CYR (I) = PYH (I) / PLR (I)
         CZR (I) = PZH (I) / PLR (I)
         ITR (I) = NREH (I)
 1000 CONTINUE
*=== End of subroutine hinhev =========================================*
      RETURN
      END
