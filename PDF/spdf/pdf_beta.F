*
* $Id$
*
* $Log$
* Revision 1.1.1.1  2000/04/19 13:29:40  mclareni
* Import version 8.04
*
* Revision 1.1  2000/03/23 11:39:41  pdflib
* beta function changed to pdf_beta function
*
* Revision 1.1.1.1  1996/04/12 15:29:04  plothow
* Version 7.01
*
*
#include "pdf/pilot.h"
C
      FUNCTION PDF_BETA(X,Y)
      PDF_BETA=GAMMA(X)*GAMMA(Y)/GAMMA(X+Y)
      RETURN
      END
