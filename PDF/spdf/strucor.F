*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:51  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:30:06  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      SUBROUTINE STRUCOR(DX,DQ,DUPV,DDNV,DSEA,DSTR,DGL)
C.......................................................................
C.    HIGH ENERGY STRUCTURE FUNCTIONS
C.
C :::::::::::: Owens and Reya   ::::::::::::::::::::::::::::
C.    Phys. Rev. D 17 (1978) 3003.
C.
C.    ( LAMBDA = 500 MEV)
C.......................................................................
#include "pdf/impdp.inc"
C
      QQ = DQ*DQ
      CALL RXUD (DX,QQ,XUD)
      CALL RXDV (DX,QQ,XDV)
      CALL RXQBAR (DX,QQ,XQBAR)
      CALL RXGLUE (DX,QQ,XGL)
      DUPV=XUD-XDV
      DDNV=XDV
      DSEA=XQBAR
      DSTR=DSEA
      DGL =XGL
      RETURN
      END
