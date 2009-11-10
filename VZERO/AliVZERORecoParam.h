#ifndef ALIVZERORECOPARAM_H
#define ALIVZERORECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with VZERO reconstruction parameters                                //
// Origin: Brigitte Cheynis                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDetectorRecoParam.h"

class AliVZERORecoParam : public AliDetectorRecoParam
{
 public: 
  AliVZERORecoParam();
  virtual ~AliVZERORecoParam();

  void SetNSigmaPed(Float_t nSigma) { fNSigmaPed = nSigma; }
  void SetStartClock(Int_t start) { fStartClock = start; }
  void SetEndClock(Int_t end) {fEndClock = end; }
  void SetNPreClocks(Int_t preClocks) { fNPreClocks = preClocks; }
  void SetNPostClocks(Int_t postClocks) { fNPostClocks = postClocks; }

  Float_t GetNSigmaPed() const { return fNSigmaPed; }
  Int_t  GetStartClock() const { return fStartClock; }
  Int_t  GetEndClock() const { return fEndClock; }
  Int_t  GetNPreClocks() const { return fNPreClocks; }
  Int_t  GetNPostClocks() const { return fNPostClocks; }

 private:

  Float_t fNSigmaPed;  // Number of pedestal sigmas for adc cut
  Int_t fStartClock;   // Start clock for max adc search
  Int_t fEndClock;     // End clock for max adc search
  Int_t fNPreClocks;   // Number of pre-clocks used in adc charge sum
  Int_t fNPostClocks;  // Number of post-clocks used in adc charge sum

  AliVZERORecoParam(const AliVZERORecoParam & param);
  AliVZERORecoParam & operator=(const AliVZERORecoParam &param);

  ClassDef(AliVZERORecoParam,2) // VZERO reco parameters
};

#endif
