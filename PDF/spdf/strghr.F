*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:31:02  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:50  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C-----------------------------------------------------------------------
      SUBROUTINE STRGHR(DX,DSCALE,DUPV,DDNV,DSEA,DSTR,DCHM,DGL)
C ::::::::::::   GHR   ::::::::::::::::::::::::::::
C
#include "pdf/expdp.inc"
     +       DX,DSCALE,DUPV,DDNV,DSEA,DSTR,DCHM,DGL
      DATA ALAM2/0.160/, Q02/4.0/
      Q = DSCALE
      X = DX
      QQ=Q*Q
      SQ= LOG( LOG( MAX(Q02,QQ)/ALAM2)/ LOG(Q02/ALAM2))
      UPV=GHRUV(X,SQ)
      DNV=GHRDV(X,SQ)
      SEA=GHRXI(X,SQ)
      STR=GHRXS(X,SQ)
      GL =GHRGL(X,SQ)
      DUPV=UPV
      DDNV=DNV
      DSEA=SEA
      DSTR=STR
      DCHM=0.0D0
      DGL =GL
      RETURN
      END
