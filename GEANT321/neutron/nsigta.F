*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:58  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE NSIGTA(E,JMED,TSIG,D,ISIGTS,LSIGT)
C       THIS ROUTINE DETERMINES THE MACROSCOPIC TOTAL
C       CROSS SECTION FOR MEDIA MED
#include "geant321/mmicab.inc"
      DIMENSION D(*),ISIGTS(*),LSIGT(*)
      CALL GTMED(JMED,MED)
      TSIG=0.0
      L1=LSIGT(MED)
      LS1=ISIGTS(MED)+LMOX3
      LEN=L1/2
      CALL TBSPLT(D(LS1),E,LEN,TSIG)
      RETURN
      END
