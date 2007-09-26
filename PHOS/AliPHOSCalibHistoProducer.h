#ifndef ALIPHOSCALIBHISTOPRODUCER_H
#define ALIPHOSCALIBHISTOPRODUCER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
// Class AliPHOSCalibHistoProducer accumulating histograms
// with amplitudes per PHOS channel
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TH1F;
class TFile;
class AliPHOSRawDecoder;

class AliPHOSCalibHistoProducer : public TObject {
public:

  AliPHOSCalibHistoProducer();
  AliPHOSCalibHistoProducer(Int_t nbinsx, Double_t xlow, Double_t xup);
  AliPHOSCalibHistoProducer(const AliPHOSCalibHistoProducer &histoproducer);
  AliPHOSCalibHistoProducer& operator= (const AliPHOSCalibHistoProducer &histoproducer);
  virtual ~AliPHOSCalibHistoProducer();

  void Run();
  void UpdateHistoFile();
  void SetUpdatingRate(Int_t rate) {fUpdatingRate = rate;}
  void SetOldRCUFormat(Bool_t isOldRCUFormat) { fIsOldRCUFormat = isOldRCUFormat; }
  void SetRawDecoder(AliPHOSRawDecoder* decoder) { fRawDecoder = decoder; }

protected:

  TH1F* fAmpHisto[5][56][64]; // amplitudes in [module][column][row].
  AliPHOSRawDecoder* fRawDecoder;   // raw data decoder.
  TFile* fHistoFile;          // root file to store histograms in
  Int_t fUpdatingRate;        // update rate
  Bool_t fIsOldRCUFormat;     // Old RCU format flag.
  Int_t fEvents;
  Int_t fNbins;               // Number of bins in histograms.
  Double_t fXlow;             // Low X in histograms.
  Double_t fXup;              // High X in histograms.

  ClassDef(AliPHOSCalibHistoProducer,1)

};

#endif
