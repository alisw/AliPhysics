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
 * Author: Jussi Viinikainen <jussi.viinikainen@cern.ch> (adaptation to Run2 format)
 */
class AliPHOSTRURawReader : public TObject
{
 public:
  AliPHOSTRURawReader();
  virtual ~AliPHOSTRURawReader();

  Short_t GetTriggerSignal(Int_t xIdx, Int_t zIdx, Int_t timeBin) const {return fSignals[xIdx][zIdx][timeBin];}
  Bool_t GetTriggerFlag(Int_t xIdx, Int_t zIdx, Int_t timeBin) const;
  bool IsActive() const {return fActive;}
  bool IsActive(Int_t timeBin) const {return fActiveTime[timeBin];}
  bool HasSignal() const {return fHasSignal;}
  bool HasSignal(Int_t timeBin) const {return fHasSignalTime[timeBin];}
  
  void ReadFromStream(AliCaloRawStreamV3* rawStream);
  void Reset();
  
  static Int_t GetDefaultSignalValue() { return fgkDefaultSignalValue; };
  
 private:
  AliPHOSTRURawReader(const AliPHOSTRURawReader &); // not implemented
  AliPHOSTRURawReader& operator= (const AliPHOSTRURawReader &); // not implemented
  
  // constants
  static const Int_t fgkDefaultSignalValue = 512; // Default/Ideal TRU amplitude pedestal
  static const Int_t fgkNReadoutChannels = 112;  // number of readout channels in tru
  static const Int_t fgkN4x4TriggerFlags = 91;  // number of possible 4x4 area in PHOS
  static const Int_t fgkFinalProductionChannel = 123; // The last channel of production bits, contains markesr to choose between 2x2 and 4x4 algorithm
  static const Int_t fgkWordLength = 10; // Length of one data word in raw data
  static const Int_t fgkNTimeBins = 128; // Number of timeBins
  static const Int_t fgkN2x2XPrTRURow = 8; // (=64/2/4) Number of 2x2 pr. row
  static const Int_t fgkN2x2ZPrBranch = 14; // (=56/2/2) Number of 2x2 pr. branch
  static const Int_t fgkN4x4XPrTRURow = 7; // (=64/2/4 -1) Number of 4x4 pr. row
  static const Int_t fgkN4x4ZPrBranch = 13; // (=56/2/2) -1 Number of 4x4 pr. branch
  
  Short_t fSignals[fgkN2x2XPrTRURow][fgkN2x2ZPrBranch][fgkNTimeBins]; // 2x2 Trigger Signal Sum, [x][z][t]
  Bool_t  fFlags[fgkN4x4XPrTRURow][fgkN4x4ZPrBranch][fgkNTimeBins]; // 4x4 Trigger Flag, [x][z][t]
  Bool_t  fFlags2x2[fgkN2x2XPrTRURow][fgkN2x2ZPrBranch][fgkNTimeBins]; // 2x2 Trigger Flag, [x][z][t]
  
  Bool_t fActive; // Active
  Bool_t fHasSignal; // Has Signal
  Bool_t fActiveTime[fgkNTimeBins]; // Active [t]
  Bool_t fHasSignalTime[fgkNTimeBins]; // Has Signal [t]
  Bool_t fUse4x4Flags; // True is the 4x4 analysis bit (bit number 112) is 1, otherwise false
  
  ClassDef(AliPHOSTRURawReader, 2)
};

#endif // ALIPHOSTRURAWREADER_H
