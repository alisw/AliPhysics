#ifndef ROOT_MevSimCOMMON
#define ROOT_MevSimCOMMON


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

enum {
  NPID = 30,            
  NFLOWTERMS = 6,       
  FACTORIAL_MAX = 10000,
  NMAXINTEG = 100,
  NMULTMAXSTEPS = 1000
};

///////////////////////////////////////////////////////////////////////////////
//
// NPID       : max # of particle ID types
// NFLOWTERMS : max # of terms in the anisotropic flow model for azimuthal
//              (phi) angle dependence. 
// FACTORIAL_MAX : max # multiplicity per event;
//                 for any specific particle ID; also used for log(n!)
// NMAXINTEG     : max # integration steps in parameter variance calculation
// NMULTMAXSTEPS : max # integration steps in multiplicity variance calculation
//                 this must be an even integer.
// Common/track/pout(npid,4,factorial_max) 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#define f2cFortran
#ifndef __CFORTRAN_LOADED
#include "cfortran.h"
#endif
#endif

#ifndef __CINT__ 
extern "C" {
  typedef struct {
//    Float_t pout[NPID][4][FACTORIAL_MAX];
     Float_t pout[NPID*4*FACTORIAL_MAX];
  }  TrackCommon;
  
#define TRACK COMMON_BLOCK(TRACK,track)
COMMON_BLOCK_DEF(TrackCommon,TRACK);
}
#endif 

#endif




