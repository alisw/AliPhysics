#ifndef ALIAMPTRNDM_H
#define ALIAMPTRNDM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <Rtypes.h>
#include <TError.h>

class TRandom;

class AliAmptRndm {
 public:
  AliAmptRndm() {}
  virtual ~AliAmptRndm() {
    fgAmptRandom=0;
  }
  
  static void SetAmptRandom(TRandom *ran=0);
  static TRandom * GetAmptRandom();

private:
  AliAmptRndm(const AliAmptRndm &Ampt);
  AliAmptRndm &operator=(const AliAmptRndm &rhs);

  static TRandom * fgAmptRandom; //! pointer to the random number generator

  ClassDef(AliAmptRndm,0)  //Random Number generator wrapper (non persistent)
};

#endif 

