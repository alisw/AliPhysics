*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#ifndef CERNLIB_GEANT321_MPSTOR_INC
#define CERNLIB_GEANT321_MPSTOR_INC
*
*
* mpstor.inc
*
      PARAMETER(IDNEU  = 1)
      PARAMETER(IDHEVY = 2)
      PARAMETER(IDGAMA = 3)
      COMMON/ MPSTOR / EN(MAXPAR),UN(MAXPAR),VN(MAXPAR),WN(MAXPAR),
     +                 AGEN(MAXPAR),MTN(MAXPAR),AMN(MAXPAR),
     +                 ZMN(MAXPAR),IDN(MAXPAR),
     +                 EP,UP,VP,WP,MTP,AGEP,AMP,ZMP,
     +                 NNEU,NHEVY,NGAMA,NPSTOR

#endif
