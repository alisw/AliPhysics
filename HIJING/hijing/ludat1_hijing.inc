#ifdef __APPLE__
#ifndef __INTEL_COMPILER
#define STOP CALL EXIT !
#define stop CALL EXIT !
#endif
#endif
* $Id$
      REAL PARU, PARJ
      COMMON/LUDAT1_HIJING/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
#ifndef BLANKET_SAVE
      SAVE /LUDAT1_HIJING/ 
#endif
