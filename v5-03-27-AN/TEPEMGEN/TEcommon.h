#ifndef ROOT_TEcommon
#define ROOT_TEcommon
//------------------------------------------------------------------------
// TEcommon is an interface COMMON blocks of the fortran event generator of
// single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 GeV/c
//%
// Yuri.Kharlov@cern.ch
// 9 October 2002
//------------------------------------------------------------------------

#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif

// c++ interface to the f77 program - event generator of
// e+e- pair production in ultraperipheral ion collisions
// Author: Yuri Kharlov, 20 September 2002

extern "C" {
  typedef struct {
    Double_t Xsect2;
    Double_t Dsect2;
    Double_t Xsecttot;
    Double_t Dsecttot;
    Int_t    Nevnt;
  } EeventCommon;
#define EEVENT COMMON_BLOCK(EEVENT,eevent)
  COMMON_BLOCK_DEF(EeventCommon,EEVENT);
}

#endif
