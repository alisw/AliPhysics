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
      FUNCTION GHRUV(Y,SQ)
C
C  STRUCTURE FUNCTION U-QUARK FROM A PROTON
      N=2
      A=0.421-0.0412*SQ
      ASQ=ABS(SQ)
      SI=1.
      IF (ASQ.GT.0.) SI=SQ/ABS(SQ)
      C=2-0.6223*SI*ASQ**0.8
      D=3.37+0.4319*SQ
      GHRUV=N*C*Y**A*(1.-Y**C)**D/PDF_BETA(D+1.0,A/C)
      RETURN
      END
