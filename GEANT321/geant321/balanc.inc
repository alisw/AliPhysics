*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:37  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_BALANC_INC
#define CERNLIB_GEANT321_BALANC_INC
*
*
* balanc.inc
*
*                                                                     *
*=== balanc ==========================================================*
*                                                                     *
      LOGICAL LRESMP, LEVDIF, LPRDIF
      COMMON /FKBALA/  ETTOT,  PTTOT, PXTTOT, PYTTOT, PZTTOT,
     &                 EINCP,  EINCN, TVGREY, TVGRE0,  TVEUZ,
     &                   EUZ,    PUX,    PUY,    PUZ,  ENUCR, PXNUCR,
     &                PYNUCR, PZNUCR,  EINTR, PXINTR, PYINTR, PZINTR,
     &                  EFRM,  PXFRM,  PYFRM,  PZFRM,   PSEA,
     &                NGREYP, NGREYN,    ICU,    IBU, ICNUCR, IBNUCR,
     &                ICINTR, IBINTR, LRESMP, LEVDIF, LPRDIF
 

#endif
