#ifndef ALIEMCALCALIBHISTOPRODUCER_H
#define ALIEMCALCALIBHISTOPRODUCER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1  2006/12/07 16:32:16  gustavo
 * First shuttle code, online calibration histograms producer, EMCAL preprocessor
 * 
 *
*/

///////////////////////////////////////////////////////////////////////////////
// Class AliEMCALCalibHistoProducer accumulating histograms
// with amplitudes per EMCAL channel
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TProfile;
class TString;
class TH1F;
class TFile;
class AliRawReader;

class AliEMCALCalibHistoProducer : public TObject {
public:

  AliEMCALCalibHistoProducer();
  AliEMCALCalibHistoProducer(AliRawReader* rawReader);
  AliEMCALCalibHistoProducer(const AliEMCALCalibHistoProducer &histoproducer);
  AliEMCALCalibHistoProducer& operator= (const AliEMCALCalibHistoProducer &histoproducer);
  virtual ~AliEMCALCalibHistoProducer();

  void Init();
  void Run();
  void UpdateHistoFile();
  void SetUpdatingRate(Int_t rate) {fUpdatingRate = rate;}
  void SetCalibHistoFileName(Int_t name) {fHistoFileName = name;}
  void SetSMInstalled(Int_t nsm, Bool_t bsm) {fSMInstalled[nsm] = bsm;}

protected:

  TH1F* fAmpHisto[12][48][24]; // amplitudes in [module][column][row].
  TProfile *fAmpProf[12]; // one per SuperModule
  AliRawReader* fRawReader;   // raw data reader.
  TFile* fHistoFile;          // root file to store histograms in
  TString fHistoFileName;          // name of root file to store histograms in
  Int_t fUpdatingRate;        // update rate
  Int_t fNSuperModules;          //Number of SuperModules;
  Int_t fNCellsEta;                  //Number of Cells in Eta in a SuperModule;
  Int_t fNCellsPhi;                   //Number of Cells in Phi in a SuperModule;
  Int_t fNCellsPhiHalfSM;      //Number of Cells in Phi in a Half SuperModule;
  Bool_t fSMInstalled[12];  //Check which detectors are on.
  ClassDef(AliEMCALCalibHistoProducer,2)

};

#endif
