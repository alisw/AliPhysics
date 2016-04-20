#ifndef ALIPHOSTRIGGERRAWDIGIT_H
#define ALIPHOSTRIGGERRAWDIGIT_H

/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliDigitNew.h"

class AliPHOSTriggerRawDigit : public AliDigitNew 
{
  
 public:

  AliPHOSTriggerRawDigit();
  AliPHOSTriggerRawDigit(Int_t module, Int_t x, Int_t z, Int_t L1_Threshold, Int_t amp);
  AliPHOSTriggerRawDigit(Int_t module, Int_t xIdx, Int_t zIdx, Int_t TRURow, Int_t branch, Int_t amp);
  AliPHOSTriggerRawDigit(const AliPHOSTriggerRawDigit & tdigit);
  AliPHOSTriggerRawDigit& operator=(const AliPHOSTriggerRawDigit & tdigit);
  virtual ~AliPHOSTriggerRawDigit() {}

  void Get4x4Position(Int_t& module, Int_t& xIdx, Int_t& zIdx, Int_t& TRURow, Int_t& branch) 
  {module = fMod; xIdx = fXIdx; zIdx = fZIdx; TRURow = fTRURow; branch = fBranch; }

  void GetModXZ(Int_t& mod, Int_t& modX, Int_t& modZ);
  Int_t GetL1Threshold() { return fL1Threshold; }
  Int_t GetType() { return fType; } 
  
 private:

  Int_t fType; // 0-L0, 1-L1
  Int_t fMod;  // module
  Int_t fXloc; // local X in module
  Int_t fZloc; // local Z in module
  Int_t fXIdx; // 4x4 X 
  Int_t fZIdx; // 4x4 Z
  Int_t fTRURow; // TRU row
  Int_t fBranch; // branch
  Int_t fL1Threshold; // L1 threshold: 0-High, 1-Medium, 2-Low
  
  ClassDef(AliPHOSTriggerRawDigit,2)
};

#endif
