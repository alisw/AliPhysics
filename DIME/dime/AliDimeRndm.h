#ifndef ALIDIMERNDM_H
#define ALIDIMERNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TRandom3.h>
#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliDimeRndm {
 public:
  AliDimeRndm() {
    // Default constructor. The static data member is initialized
    // in the implementation file
  }
  AliDimeRndm(const AliDimeRndm &/*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliDimeRndm() {
    // Destructor
    fgDimeRandom=0;
  }
  AliDimeRndm & operator=(const AliDimeRndm& /*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetDimeRandom(TRandom3* ran = 0);
  static TRandom3* GetDimeRandom();

private:

  static TRandom3* fgDimeRandom; //! pointer to the random number generator

  ClassDef(AliDimeRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

