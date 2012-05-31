#ifndef ALIPHOSTRIGGERRAWDIGIT_H
#define ALIPHOSTRIGGERRAWDIGIT_H

/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

#include "AliDigitNew.h"

class AliPHOSTriggerRawDigit : public AliDigitNew 
{
 public:

  AliPHOSTriggerRawDigit();
  AliPHOSTriggerRawDigit(Int_t module, Int_t xIdx, Int_t zIdx, Int_t TRURow, Int_t branch, Int_t amp);
  AliPHOSTriggerRawDigit(const AliPHOSTriggerRawDigit & tdigit);
  AliPHOSTriggerRawDigit& operator=(const AliPHOSTriggerRawDigit & tdigit);
  virtual ~AliPHOSTriggerRawDigit() {}

  void Get4x4Position(Int_t& module, Int_t& xIdx, Int_t& zIdx, Int_t& TRURow, Int_t& branch) 
  {module = fMod; xIdx = fXIdx; zIdx = fZIdx; TRURow = fTRURow; branch = fBranch; }

  void GetModXZ(Int_t& mod, Int_t& modX, Int_t& modZ);

 private:

  Int_t fMod;  // module
  Int_t fXIdx; // 4x4 X 
  Int_t fZIdx; // 4x4 Z
  Int_t fTRURow; // TRU row
  Int_t fBranch; // branch

  ClassDef(AliPHOSTriggerRawDigit,1)
};

#endif
