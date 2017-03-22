#ifndef ALIPYTHIARNDM_H
#define ALIPYTHIARNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliPythiaRndm.h 8920 2004-01-13 11:29:51Z hristov $ */

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliPythiaRndm {
 public:
  AliPythiaRndm() {
    // Default constructor. The static data member is initialized 
    // in the implementation file
  }
  AliPythiaRndm(const AliPythiaRndm & /*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliPythiaRndm() {
    // Destructor
    fgPythiaRandom=0;
  }
  AliPythiaRndm & operator=(const AliPythiaRndm& /*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetPythiaRandom(TRandom *ran=0);
  static TRandom * GetPythiaRandom();

private:

  static TRandom * fgPythiaRandom; //! pointer to the random number generator

  ClassDef(AliPythiaRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

