*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:36  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE GSCANO
*
************************************************************************
*             SCAN utility called at the end of each event.            *
*                                                                      *
*          ==>Called by : <USER>,GUKINE                                *
*             Author    R.Brun  *********                              *
************************************************************************
*                                                                      *
#include "geant321/gcbank.inc"
#include "geant321/gcscal.inc"
#include "geant321/gcscan.inc"
C.
C.    ------------------------------------------------------------------
C.
*
*          Insert the temporary LSTEMP structure into LSCAN
*
      IF(IPHI.GT.IQ(LSCAN-2))GO TO 10
      IF(LQ(LSCAN-IPHI).NE.0)GO TO 10
      NPUSH=ISCUR-IQ(LSTEMP-1)-1
      IQ(LSTEMP-5)=IPHI
      CALL MZPUSH(IXCONS,LSTEMP,0,NPUSH,'I')
      CALL ZSHUNT(IXCONS,LSTEMP,LSCAN,-IPHI,0)
      GO TO 999
   10 IF(LSTEMP.NE.0)CALL MZDROP(IXCONS,LSTEMP,' ')
      LSTEMP=0
*
  999 END
