#ifndef ALIDPMJETRNDM_H
#define ALIDPMJETRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliDpmJetRndm {
 public:
  AliDpmJetRndm() {
    // Default constructor. The static data member is initialized 
    // in the implementation file
  }
  AliDpmJetRndm(const AliDpmJetRndm & /*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliDpmJetRndm() {
  // Destructor
    fgDpmJetRandom=0;
  }

  AliDpmJetRndm & operator=(const AliDpmJetRndm& /*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetDpmJetRandom(TRandom *ran=0);
  static TRandom * GetDpmJetRandom();

private:

  static TRandom * fgDpmJetRandom; //! pointer to the random number generator

  ClassDef(AliDpmJetRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

