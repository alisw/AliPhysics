#ifndef ALIHIJINGRNDM_H
#define ALIHIJINGRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliHijingRndm {
 public:
  AliHijingRndm() {
    // Default constructor. The static data member is initialized
    // in the implementation file
  }
  AliHijingRndm(const AliHijingRndm &/*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliHijingRndm() {
    // Destructor
    fgHijingRandom=0;
  }
  AliHijingRndm & operator=(const AliHijingRndm& /*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetHijingRandom(TRandom *ran=0);
  static TRandom * GetHijingRandom();

private:

  static TRandom * fgHijingRandom; //! pointer to the random number generator

  ClassDef(AliHijingRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

