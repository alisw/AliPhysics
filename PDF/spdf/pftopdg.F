*
* $Id$
*
* $Log$
* Revision 1.1.1.2  1996/10/30 08:29:06  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:26  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      SUBROUTINE PFTOPDG(DX,DSCALE,DXPDF)
C
#include "pdf/expdp.inc"
     +       DX,DSCALE,DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL,
     +       DXPDF(-6:6)
C... call STRUCTM in PDFLIB to get flavour content
      CALL STRUCTM(DX,DSCALE,
     +                   DUPV,DDNV,DUSEA,DDSEA,DSTR,DCHM,DBOT,DTOP,DGL)
C... convert flavour convention of PDFLIB to PDG convention
      DXPDF(0) = DGL
      DXPDF(1) = DDNV + DDSEA
      DXPDF(2) = DUPV + DUSEA
      DXPDF(3) = DSTR
      DXPDF(4) = DCHM
      DXPDF(5) = DBOT
      DXPDF(6) = DTOP
      DXPDF(-1) = DDSEA
      DXPDF(-2) = DUSEA
      DXPDF(-3) = DSTR
      DXPDF(-4) = DCHM
      DXPDF(-5) = DBOT
      DXPDF(-6) = DTOP
C
      RETURN
      END
