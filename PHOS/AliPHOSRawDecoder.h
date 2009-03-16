#ifndef ALIPHOSRAWDECODER_H
#define ALIPHOSRAWDECODER_H
/* Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                          */

/* $Id$ */

// This class extracts the PHOS "digits" of current event
// (amplitude,time, position,gain) from the raw stream 
// provided by AliRawReader. See cxx source for use case.

#include "AliCaloRawStream.h"

class TArrayI;
class AliRawReader;
class AliPHOSCalibData ;

class AliPHOSRawDecoder : public TObject 
{

public:

  AliPHOSRawDecoder();
  AliPHOSRawDecoder(AliRawReader* rawReader, AliAltroMapping **mapping = NULL);
  AliPHOSRawDecoder(const AliPHOSRawDecoder& rawDecoder);
  AliPHOSRawDecoder& operator = (const AliPHOSRawDecoder& rawDecoder);
  virtual ~AliPHOSRawDecoder();

  virtual Bool_t NextDigit();

  void SubtractPedestals(Bool_t subtract) {fPedSubtract=subtract;}
  void SetAmpOffset(Int_t extPed=5){fAmpOffset=extPed ;}
  void SetAmpThreshold(Int_t thr=5){fAmpThreshold=thr ;}

  Double_t GetEnergy() const { return fEnergy; }
  Double_t GetTime() const { return fTime; }
  Double_t GetSampleQuality() const {return fQuality ;}
  Double_t GetPedestalRMS() const {return fPedestalRMS ;}
  Int_t GetModule() const { return fModule; }
  Int_t GetColumn() const { return fColumn; }
  Int_t GetRow() const { return fRow; }
  Bool_t IsLowGain() const { return fLowGainFlag; }
  Bool_t IsOverflow() const { return fOverflow ;}

  const AliRawReader* GetRawReader() const { return fRawReader; }
  void SetCalibData(AliPHOSCalibData * cdata){ fCalibData=cdata ;}

protected:   
  
  AliRawReader* fRawReader;      // raw data reader
  AliCaloRawStream* fCaloStream; // PHOS/EMCAL raw data stream
  Bool_t fPedSubtract;           // pedestals subtraction (kTRUE="yes")


  Double_t fEnergy; // "digit" energy
  Double_t fTime;   // "digit" time
  Double_t fQuality ; //Sample quality
  Double_t fPedestalRMS; //calciulated RMS of pedestal (non-ZS runs)
  Int_t fAmpOffset ; //Pedestal offset from ALTRO chips
  Int_t fAmpThreshold ; //Zero Suppression threshold from ALTRO chips
  Int_t fModule;    // PHOS module number (1-5)
  Int_t fColumn;    // column in the module
  Int_t fRow;       // row
  Int_t fNewModule;    // PHOS module number (1-5) of keeped sample
  Int_t fNewColumn;    // column in the module  of keeped sample
  Int_t fNewRow;       // row  of keeped sample
  Int_t fNewAmp ;      //Keeped amp
  Int_t fNewTime ;     //Time of keeped sample
  Bool_t fLowGainFlag; //True if sample read from Low Gain
  Bool_t fNewLowGainFlag; // fLowGainFlag of keeped sample
  Bool_t fOverflow ;   //Wether there was overflow
  TArrayI* fSamples;   // array of samples
  TArrayI* fTimes ;    // array of times corresponding to samples
  AliPHOSCalibData * fCalibData ;   //! Calibration database if avalable


  ClassDef(AliPHOSRawDecoder,4)
};

#endif
