#ifndef ALIPHOSRAWDECODER_H
#define ALIPHOSRAWDECODER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliRawReader.h"
#include "AliCaloRawStream.h"

class AliPHOSRawDecoder {

public:

  AliPHOSRawDecoder();
  AliPHOSRawDecoder(AliRawReader* rawReader);
  AliPHOSRawDecoder(const AliPHOSRawDecoder& rawDecoder);
  AliPHOSRawDecoder& operator = (const AliPHOSRawDecoder& rawDecoder);
  virtual ~AliPHOSRawDecoder();

  virtual Bool_t NextDigit();

  void SetOldRCUFormat(Bool_t isOldRCU) {fCaloStream->SetOldRCUFormat(isOldRCU);}
  void SubtractPedestals(Bool_t subtract) {fPedSubtract=subtract;}

  Double_t GetEnergy() { return fEnergy; }
  Double_t GetTime() { return fTime; }
  Int_t GetModule() { return fModule; }
  Int_t GetColumn() { return fColumn; }
  Int_t GetRow() { return fRow; }
  Bool_t IsLowGain() { return fCaloStream->IsLowGain(); }

protected:   
  
  AliRawReader* fRawReader;
  AliCaloRawStream* fCaloStream;
  Bool_t fPedSubtract;

private:

  Double_t fEnergy;
  Double_t fTime;
  Int_t fModule;
  Int_t fColumn;
  Int_t fRow;
  
  ClassDef(AliPHOSRawDecoder,1)
};

#endif
