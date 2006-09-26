#ifndef ALIISAJETRNDM_H
#define ALIISAJETRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliIsajetRndm {
 public:
  AliIsajetRndm() {
    // Default constructor. The static data member is initialized
    // in the implementation file
  }
  AliIsajetRndm(const AliIsajetRndm &/*rn*/) {
    // Copy constructor: no copy allowed for the object
    ::Fatal("Copy constructor","Not allowed\n");
  }
  virtual ~AliIsajetRndm() {
    // Destructor
    fgIsajetRandom=0;
  }
  AliIsajetRndm & operator=(const AliIsajetRndm &/*rn*/) {
    // Assignment operator: no assignment allowed
    ::Fatal("Assignment operator","Not allowed\n");
    return (*this);
  }
  
  static void SetIsajetRandom(TRandom *ran=0);
  static TRandom * GetIsajetRandom();

private:

  static TRandom * fgIsajetRandom; //! pointer to the random number generator

  ClassDef(AliIsajetRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

