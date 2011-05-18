#ifndef AliTRDPTRGCBB_H
#define AliTRDPTRGCBB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// --------------------------------------------------------
// 
// PTRG simulation
//
// --------------------------------------------------------
#include "AliTRDptrgParam.h"
#include <TObjArray.h>
#include <TObject.h>

class AliRunLoader;

class AliTRDptrgCBAC;
class AliTRDptrgTLMU;

class AliTRDptrgCBB : public TObject {
 public:
  AliTRDptrgCBB(AliRunLoader *rl = 0x0);
  AliTRDptrgCBB(AliRunLoader *rl, AliTRDptrgParam* param, 
                AliTRDptrgOperatingMode_t operatingMode);
  ~AliTRDptrgCBB();
  
  Int_t* Simulate(); // Simulates the ptrg behavior of event
  Bool_t GetPT(); // Evaluates ptrg decision
 protected:
  Bool_t LoadParams(); // loads the parameters stored

  AliRunLoader *fRunLoader;  //!
  AliTRDptrgParam *fParam; // singleton obj containing configuration parameters
  AliTRDptrgOperatingMode_t fOperatingMode; // working on Digits or Hits?

  AliTRDptrgCBAC *fCBA; // control box at a side of the solenoid
  AliTRDptrgCBAC *fCBC; // control box at c side of the solenoid
  AliTRDptrgTLMU *fTLMU; // TLMU

  TObjArray fLUTArray; // Array with Look-Up-Tables (usually two, called X,Y)

  const AliTRDptrgParam::AliTRDptrgPTmasks *fPTmasks; // PT output masks 
 private:
  AliTRDptrgCBB& operator=(const AliTRDptrgCBB &rhs); // not implemented
  AliTRDptrgCBB(const AliTRDptrgCBB &rhs); // not implemented

  ClassDef(AliTRDptrgCBB, 1);
};

#endif
