*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:19:35  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.45  by  S.Giani
*-- Author :
*$ CREATE BLKDT2.FOR
*COPY BLKDT2
*
*=== blkdt2 ===========================================================*
*==                                                                    *
      SUBROUTINE FLKDT2
 
#include "geant321/dblprc.inc"
#include "geant321/dimpar.inc"
#include "geant321/iounit.inc"
*         CHARACTER*8 EDI, EDI1, EDI2      Not used here! A.F.
*
      DIMENSION RUNT1(2)
*
      COMMON /FKCD97/ CDDT
        DATA CDDT1 / .3D0 /
      COMMON /FKRUN/ RUNT(2)
        DATA RUNT1 / 2*100.D0 /
      COMMON/FKIDGB/ IDGB, IGB
        DATA IDGB1,IGB1 / 0, 0 /
      COMMON / FKENCO/ ETEST, TNKTE
      COMMON / FKNUCT / ETHR, PTHR
        DATA ETHR1,TNKTE1 / 2*1.D-3 /
*
      CDDT    = CDDT1
      RUNT(1) = RUNT1(1)
      RUNT(2) = RUNT1(2)
      IDGB    = IDGB1
      IGB     = IGB1
      ETHR    = ETHR1
      TNKTE   = TNKTE1
*=                                               end*block.blkdt2      *
      END
