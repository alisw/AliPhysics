*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:56  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.48  by  S.Giani
*-- Author :
      SUBROUTINE GETPAR(ID,N,IERR)
C retrieve particle from MPSTOR common
#include "geant321/mconst.inc"
#include "geant321/mpstor.inc"
#include "geant321/minput.inc"
      IERR = 0
      NN = 0
      NS = 1
   10 CONTINUE
      IF(IDN(NS).EQ.ID) NN = NN + 1
      IF(N.EQ.NN) GOTO 20
      NS = NS + 1
      IF(NS.GT.NPSTOR) THEN
         WRITE(IOUT,'('' MICAP: Cant retrieve particle no. '',I3,      '
     +   //'          '' of type '',I3,''; End of data '')') N,ID
         IERR = 1
         RETURN
      ENDIF
      GOTO 10
   20 CONTINUE
      EP = EN(NS)
      UP = UN(NS)
      VP = VN(NS)
      WP = WN(NS)
      AMP = AMN(NS)
      ZMP = ZMN(NS)
      AGEP = AGEN(NS)
      MTP = MTN(NS)
      RETURN
      END
