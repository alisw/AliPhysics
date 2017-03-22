#ifndef ALIHERWIGRNDM_H
#define ALIHERWIGRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------------
//   Class: AliHerwigRndm
//   Responsibilities: Interface to Root random number generator 
//                     from Fortran (re-implements FINCTION RLU_HERWIG 
//                     from HERWIG)
//-----------------------------------------------------------------------------

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliHerwigRndm {
 public:
  AliHerwigRndm() {
    // Default constructor. The static data member is initialized
    // in the implementation file
  }
  AliHerwigRndm(const AliHerwigRndm &/*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliHerwigRndm() {
    // Destructor
    fgHerwigRandom=0;
  }
  AliHerwigRndm & operator=(const AliHerwigRndm& /*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetHerwigRandom(TRandom *ran=0);
  static TRandom * GetHerwigRandom();

private:

  static TRandom * fgHerwigRandom; //! pointer to the random number generator

  ClassDef(AliHerwigRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

