*
* $Id$
*
* $Log$
* Revision 1.1.1.3  2000/04/19 13:28:08  mclareni
* Import version 8.04
*
* Revision 1.2  2000/03/23 11:39:41  pdflib
* beta function changed to pdf_beta function
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
      FUNCTION GHRDV(Y,SQ)
C
C  D-QUARK FROM A PROTON STRUCTURE FUNCTION
      N=1
      A=0.364-0.0368*SQ
      ASQ=ABS(SQ)
      SI=1.
      IF (ASQ.GT.0.) SI=SQ/ABS(SQ)
      C=2-0.5414*SI*ASQ**0.8
      D=5.09+0.3463*SQ
      GHRDV=N*C*Y**A*(1.-Y**C)**D/PDF_BETA(D+1.0,A/C)
      RETURN
      END
