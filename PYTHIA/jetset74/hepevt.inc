*
* $Id$
*
* $Log$
* Revision 1.2  1996/05/24 16:02:35  cernlib
* Add Double Precision  PHEP,VHEP
*
* Revision 1.1.1.1  1996/03/08 17:32:19  mclareni
* jetset74
*
*
#ifndef CERNLIB_JETSET74_HEPEVT_INC
#define CERNLIB_JETSET74_HEPEVT_INC
*
*
* hepevt.inc
*
#include "jetset74/nmxhep.inc"
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/

#endif
