*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:21:49  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.33  by  S.Giani
*-- Author :
      SUBROUTINE GXCLOS(LUN)
*
*                To close unit LUN
*
#if defined(CERNLIB_IBM)
      CHARACTER*80 CHTEMP
      CHARACTER*3 CHLUN
#endif
*________________________________________________________
*
         CLOSE (LUN)
#if defined(CERNLIB_IBM)
         WRITE(CHLUN,'(I3)')LUN
         CHTEMP='FILEDEF '//CHLUN//'  CLEAR'
         CALL VMCMS(CHTEMP,ISTAT)
#endif
*
      END
 
