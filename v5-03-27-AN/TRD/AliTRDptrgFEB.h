#ifndef ALITRDPTRGFEB_H
#define ALITRDPTRGFEB_H
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

class AliTRDptrgFEB : public TObject {
 public:
  AliTRDptrgFEB(AliRunLoader *rl = 0x0);
  AliTRDptrgFEB(AliRunLoader *rl, AliTRDptrgFEBType_t febType, 
                AliTRDptrgOperatingMode_t operatingMode, 
                AliTRDptrgFEBPosition_t position, Int_t id,
                AliTRDptrgParam *param);

  ~AliTRDptrgFEB();
  Int_t* Simulate(); // starts a simulation
protected:
  Int_t LoadDigits(); // loads Digits (for usage with aquired data)
  Int_t LoadAndProcessHits(); 
  // load and process hits (for usage with simulated data)
  Bool_t LoadParams(); // load AliTRDprtgParam content
  

  AliRunLoader *fRunLoader;  //!
  AliTRDptrgParam *fParam; // Configuration parameter object
  TObjArray fLUTArray; // Array with Look-Up-Tables

  AliTRDptrgFEBType_t fType; // Indicates what input FEB uses (V0 or T0)
  AliTRDptrgOperatingMode_t fOperatingMode; // working on Digits or Hits?
  Int_t fInputChannelCount; // Number of input channels 
  AliTRDptrgFEBPosition_t fPosition; // 0 = unkown, 1 = A, 2 = C
  Int_t fID; // 0 = T0, 1 = V0-1, 2 = V0-2, 3 = V0-3, 4 = V0-4 (numbering?)
  
  UInt_t *fThreshold; // specifies the threshold for incoming analog signals
 private:
  AliTRDptrgFEB& operator=(const AliTRDptrgFEB &rhs); // not implemented
  AliTRDptrgFEB(const AliTRDptrgFEB &rhs); // not implemented

  ClassDef(AliTRDptrgFEB, 1);
};

#endif
