*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:30:58  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:49  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      SUBROUTINE STRBEP(DX,DQ,DUPV,DDNV,DSEA,DSTR,DGL)
C.......................................................................
C.
C.    HIGH ENERGY STRUCTURE FUNCTIONS
C.
C :::::::::::: Baier et al.     ::::::::::::::::::::::::::::
C.    Z. Phy. C 2 (1979) 265.
C.
C.    ( LAMBDA = 470 MEV)
C.......................................................................
#include "pdf/impdp.inc"
C
      QQ = DQ*DQ
      CALL BXF3 (DX,QQ,XF3)
      CALL BXQBAR (DX,QQ,XQBAR)
      CALL BXGLUE (DX,QQ,XGL)
      DUPV=2.D0/3.D0*XF3
      DDNV=DUPV/2.D0
      DSEA=XQBAR/6.D0
      DSTR=DSEA
      DGL =XGL
      RETURN
      END
