*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:59  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE STOPAR(ID,NP)
C store particle in MPSTOR common
#include "geant321/mconst.inc"
#include "geant321/mpstor.inc"
#include "geant321/minput.inc"
C
      NPSTOR = NPSTOR + 1
      IF(NPSTOR.GT.MAXPAR) THEN
         WRITE(IOUT,'('' MICAP :  Cant store particle; bank full'',    '
     +   //'                 '' ID='',I3,'' NPSTOR='',I5)') ID,NPSTOR
         NPSTOR = NPSTOR - 1
      ELSE
         EN(NPSTOR) = EP
         UN(NPSTOR) = UP
         VN(NPSTOR) = VP
         WN(NPSTOR) = WP
         AMN(NPSTOR) = AMP
         ZMN(NPSTOR) = ZMP
         AGEN(NPSTOR) = AGEP
         MTN(NPSTOR) = MTP
         IDN(NPSTOR) = ID
         NP = NP + 1
      ENDIF
      RETURN
      END
