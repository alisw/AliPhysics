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
  void SetAdcThresHold(Float_t val) { fAdcThresHold = val;}
  void SetTimeWindowBBALow(Float_t val) { fTimeWindowBBALow = val; }
  void SetTimeWindowBBAUp (Float_t val) { fTimeWindowBBAUp  = val; }
  void SetTimeWindowBGALow(Float_t val) { fTimeWindowBGALow = val; }
  void SetTimeWindowBGAUp (Float_t val) { fTimeWindowBGAUp  = val; }
  void SetTimeWindowBBCLow(Float_t val) { fTimeWindowBBCLow = val; }
  void SetTimeWindowBBCUp (Float_t val) { fTimeWindowBBCUp  = val; }
  void SetTimeWindowBGCLow(Float_t val) { fTimeWindowBGCLow = val; }
  void SetTimeWindowBGCUp (Float_t val) { fTimeWindowBGCUp  = val; }
  void SetMaxResid (Float_t val) { fMaxResid  = val; }

  Float_t GetNSigmaPed() const { return fNSigmaPed; }
  Int_t  GetStartClock() const { return fStartClock; }
  Int_t  GetEndClock() const { return fEndClock; }
  Int_t  GetNPreClocks() const { return fNPreClocks; }
  Int_t  GetNPostClocks() const { return fNPostClocks; }
  Float_t  GetAdcThresHold() const { return fAdcThresHold; }
  Float_t  GetTimeWindowBBALow() const { return fTimeWindowBBALow; }
  Float_t  GetTimeWindowBBAUp () const { return fTimeWindowBBAUp ; }
  Float_t  GetTimeWindowBGALow() const { return fTimeWindowBGALow; }
  Float_t  GetTimeWindowBGAUp () const { return fTimeWindowBGAUp ; }
  Float_t  GetTimeWindowBBCLow() const { return fTimeWindowBBCLow; }
  Float_t  GetTimeWindowBBCUp () const { return fTimeWindowBBCUp ; }
  Float_t  GetTimeWindowBGCLow() const { return fTimeWindowBGCLow; }
  Float_t  GetTimeWindowBGCUp () const { return fTimeWindowBGCUp ; }
  Float_t  GetMaxResid () const { return fMaxResid; }

 private:

  Float_t fNSigmaPed;  // Number of pedestal sigmas for adc cut
  Int_t fStartClock;   // Start clock for max adc search
  Int_t fEndClock;     // End clock for max adc search
  Int_t fNPreClocks;   // Number of pre-clocks used in adc charge sum
  Int_t fNPostClocks;  // Number of post-clocks used in adc charge sum
  // Cuts used in the trigger mask creation
  Float_t fAdcThresHold; // Threshold on the ADC
  Float_t fTimeWindowBBALow;  // BBA window (lower cut)
  Float_t fTimeWindowBBAUp;   // BBA window (upper cut)
  Float_t fTimeWindowBGALow;  // BGA window (lower cut)
  Float_t fTimeWindowBGAUp;   // BGA window (upper cut)
  Float_t fTimeWindowBBCLow;  // BBC window (lower cut)
  Float_t fTimeWindowBBCUp;   // BBC window (upper cut)
  Float_t fTimeWindowBGCLow;  // BGC window (lower cut)
  Float_t fTimeWindowBGCUp;   // BGC window (upper cut)
  Float_t fMaxResid;   	      // Maximum residual of a single channel time

  AliVZERORecoParam(const AliVZERORecoParam & param);
  AliVZERORecoParam & operator=(const AliVZERORecoParam &param);

  ClassDef(AliVZERORecoParam,3) // VZERO reco parameters
};

#endif
