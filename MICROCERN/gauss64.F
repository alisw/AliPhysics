*
* $Id$
*
* $Log$
* Revision 1.2.4.1  2002/11/26 16:46:57  hristov
* Merging NewIO with v3-09-04
*
* Revision 1.2  2002/10/14 14:57:10  hristov
* Merging the VirtualMC branch to the main development branch (HEAD)
*
* Revision 1.1.2.1  2002/08/29 13:47:41  alibrary
* Adding missing routines
*
* Revision 1.2  2002/07/26 15:33:28  hristov
* Dummy subroutines needed by Sun CC5.2
*
* Revision 1.1.1.1  1999/05/18 15:55:34  fca
* AliRoot sources
*
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
#include "gen/pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION DGAUSS(F,A,B,EPS)
#include "gen/imp64.inc"

      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSS')
#include "gausscod.inc"

   99 DGAUSS=H
      RETURN
      END
#else
      SUBROUTINE DGAUSS_DUMMY
      END
#endif
