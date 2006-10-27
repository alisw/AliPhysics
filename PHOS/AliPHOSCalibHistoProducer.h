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
class AliRawReader;

class AliPHOSCalibHistoProducer : public TObject {
public:

  AliPHOSCalibHistoProducer();
  AliPHOSCalibHistoProducer(AliRawReader* rawReader);
  AliPHOSCalibHistoProducer(const AliPHOSCalibHistoProducer &histoproducer);
  AliPHOSCalibHistoProducer& operator= (const AliPHOSCalibHistoProducer &histoproducer);
  virtual ~AliPHOSCalibHistoProducer();

  void Run();
  void UpdateHistoFile();
  void SetUpdatingRate(const Int_t rate) {fUpdatingRate = rate;}
  void SetOldRCUFormat(Bool_t isOldRCUFormat) { fIsOldRCUFormat = isOldRCUFormat; }

protected:

  TH1F* fAmpHisto[5][56][64]; // amplitudes in [module][column][row].
  AliRawReader* fRawReader;   // raw data reader.
  TFile* fHistoFile;          // root file to store histograms in
  Int_t fUpdatingRate;        // update rate
  Bool_t fIsOldRCUFormat;     // Old RCU format flag.

  ClassDef(AliPHOSCalibHistoProducer,1)

};

#endif
