#ifndef AliTRDPTRGTLMU_H
#define AliTRDPTRGTLMU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --------------------------------------------------------
// 
// PTRG simulation
//
// --------------------------------------------------------

#include "TObject.h"
#include "AliTRDptrgParam.h"

class AliRunLoader;

class AliTRDptrgParam;

class AliTRDptrgTLMU : public TObject {
 public:
  AliTRDptrgTLMU(AliRunLoader *rl = 0x0);
  AliTRDptrgTLMU(AliRunLoader *rl,  AliTRDptrgParam *param, 
                 AliTRDptrgOperatingMode_t operatingMode);
  ~AliTRDptrgTLMU();

  Int_t* Simulate(); // starts a simulation
  

 protected:
  Bool_t LoadParams(); // load AliTRDprtgParam content
  
  // functions for input data processing ---------------------------------------
  Int_t LoadDigits(); // loads Digits (for usage with aquired data)
  void GetInputBits(); // Gets TOF-to-TRD bits from AliTOFTrigger

  // logical functions ---------------------------------------------------------
  Int_t BackToBack(Int_t iSM, Int_t range = 0); // Back-To-Back check
  // (for +-1 and so on) SM0 and SM8, SM9, SM10 (range == 1)
  Int_t Coincidence(Int_t iSM1, Int_t iSM2); // more flexible version of 
  // BackToBack(..)
  inline Int_t Or(Int_t iSM); // returns >=1 for iSM
  Int_t GetMultiplicity(Int_t iSM);  // returns multiplicity of supermodule iSM 
  Int_t GetMultiplicitySum(); // returns the multiplicity of the whole detector

  UInt_t GetBitVectorMultiplicity(UInt_t BitVector); 
  // returns the multiplicity of a bit vector

  // variables -----------------------------------------------------------------
  AliRunLoader *fRunLoader;  //!
  AliTRDptrgParam *fParam; // Configuration parameter object
  AliTRDptrgOperatingMode_t fOperatingMode; // working on Digits or Hits?
  
  const UInt_t* fInputMask; // input mask for TOF-bits (18x32=576)
  UInt_t fTOFinputBits[18]; // input bits from TOF (18x32)
  
  UInt_t** fCMatrices;    // get coincidence matrices
  UInt_t** fMultiplicity;    // get multiplicity slices
  Int_t** fOutput;    // get output signal assignment
  
 private:
  AliTRDptrgTLMU& operator=(const AliTRDptrgTLMU &rhs); // not implemented
  AliTRDptrgTLMU(const AliTRDptrgTLMU &rhs); // not implemented

  enum{
    kNLTM = 72,          //Number of LTM
    kNLTMchannels = 48,  //Number of channels in a LTM
    kNCTTM = 36,         //Number of CTTM per TOF side
    kNCTTMchannels = 24  //Number of channels in a CTTM
  };

  ClassDef(AliTRDptrgTLMU, 1);
};

#endif
