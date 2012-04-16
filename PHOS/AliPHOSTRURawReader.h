#ifndef ALIPHOSTRURAWREADER_H
#define ALIPHOSTRURAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */


class AliCaloRawStreamV3;

#include "TObject.h"

/* 
 * Class for reading TRU data from a bunch from a raw datastream.
 * Author: Henrik Qvigstad <henrik.qvigstad@cern.ch>
 */
class AliPHOSTRURawReader : public TObject
{
 public:
  AliPHOSTRURawReader();
  virtual ~AliPHOSTRURawReader();

  Short_t GetTriggerSignal(Int_t xIdx, Int_t zIdx, Int_t timeBin) const {return fSignals[xIdx][zIdx][timeBin];}
  Bool_t GetTriggerFlag(Int_t xIdx, Int_t zIdx, Int_t timeBin) const {return fFlags[xIdx][zIdx][timeBin];}
  bool IsActive() {return fActive;}
  bool IsActive(Int_t timeBin) {return fActiveTime[timeBin];}
  bool HasSignal() {return fHasSignal;}
  bool HasSignal(Int_t timeBin) {return fHasSignalTime[timeBin];}
  
  void ReadFromStream(AliCaloRawStreamV3* );
  void Reset();

 private:
  AliPHOSTRURawReader(const AliPHOSTRURawReader &); // not implemented
  AliPHOSTRURawReader& operator= (const AliPHOSTRURawReader &); // not implemented
  
 public:
  const static Int_t kDefaultSignalValue = 512; // Default/Ideal TRU amplitude pedestal
  
 private:  
  // constants
  const static Int_t kNTimeBins = 128; // Number of timeBins
  const static Int_t kN2x2XPrTRURow = 8; // (=64/2/4) Number of 2x2 pr. row
  const static Int_t kN2x2ZPrBranch = 14; // (=56/2/2) Number of 2x2 pr. branch
  const static Int_t kN4x4XPrTRURow = 7; // (=64/2/4 -1) Number of 4x4 pr. row
  const static Int_t kN4x4ZPrBranch = 13; // (=56/2/2) -1 Number of 4x4 pr. branch
  
  Short_t fSignals[kN2x2XPrTRURow][kN2x2ZPrBranch][kNTimeBins]; // 2x2 Trigger Signal Sum, [x][z][t]
  Bool_t  fFlags[kN4x4XPrTRURow][kN4x4ZPrBranch][kNTimeBins]; // 4x4 Trigger Flag, [x][z][t]
  
  Bool_t fActive; // Active
  Bool_t fHasSignal; // Has Signal
  Bool_t fActiveTime[kNTimeBins]; // Active [t]
  Bool_t fHasSignalTime[kNTimeBins]; // Has Signal [t]
  
  ClassDef(AliPHOSTRURawReader, 1)
};

#endif // ALIPHOSTRURAWREADER_H
