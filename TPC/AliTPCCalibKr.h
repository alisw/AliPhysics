#ifndef ALITPCCALIBKR_H
#define ALITPCCALIBKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliTPCclusterKr.h"

class TH3F;
class TH1D;

class AliTPCCalibKr : public TObject {

public:
  AliTPCCalibKr();
  AliTPCCalibKr(const AliTPCCalibKr&); // copy constructor
  virtual ~AliTPCCalibKr();

  AliTPCCalibKr& operator=(const AliTPCCalibKr&); 

  //
  void Init();
  Bool_t Process(AliTPCclusterKr *cluster);
  Bool_t Accept(AliTPCclusterKr *cluster);
  Bool_t Update(AliTPCclusterKr *cluster);
  TH3F*  CreateHisto(Int_t chamber);

  const TObjArray* GetHistoKrArray () {return &fHistoKrArray;}  // get calibration object
  TH3F* GetHistoKr(Int_t sector) const;                         // get refernce histogram

  Bool_t IsCSide(Int_t chamber);
  Bool_t IsIROC(Int_t chamber);

  void Analyse();
  static TH1D* ProjectHisto(TH3F* histo3D, const char* name = "_pz", Int_t firstxbin = 0, Int_t lastxbin = 0, Int_t firstybin = 0, Int_t lastybin = 0);

  void SetASide(Bool_t bA = kTRUE) {fASide = bA;} // fill histo only A TPC side
  void SetBSide(Bool_t bC = kTRUE) {fCSide = bC;} // fill histo only C TPC side

  //Merge output objects (needed by PROOF)
  virtual Long64_t Merge(TCollection* list);

  // Setters
  void SetADCOverClustSizeRange(Float_t min=0.0,Float_t max=1.0e9)   {fADCOverClustSizeMin = min ; fADCOverClustSizeMax = max;  }
  void SetMaxADCOverClustADCRange(Float_t min=0.0,Float_t max=1.0e9) {fMaxADCOverClustADCMin = min ; fMaxADCOverClustADCMax = max;  }
  void SetTimeRange(Float_t min=0.0, Float_t max=1.0e9)              {fTimeMin = min ; fTimeMax = max; }
  void SetClustSizeRange(Float_t min=0.0, Float_t max=1.0e9)         {fClustSizeMin = min ; fClustSizeMax = max; }

private:

  Bool_t fASide;              //! Only A side
  Bool_t fCSide;              //! Only C side 
  TObjArray fHistoKrArray;    //  Calibration histograms for Kr distribution

  Float_t fADCOverClustSizeMin; // min ADCcluster over Cluster size ratio
  Float_t fADCOverClustSizeMax; // max ADCcluster over Cluster size ratio
  Float_t fMaxADCOverClustADCMin; // min MaxADC over ADCcluster ratio 
  Float_t fMaxADCOverClustADCMax; // max MaxADC over ADCcluster ratio
  Float_t fTimeMin; // min time bin for MaxADC
  Float_t fTimeMax; // max time bin for MaxADC
  Float_t fClustSizeMin; // min cluster size
  Float_t fClustSizeMax; // max cluster size

public:
  ClassDef(AliTPCCalibKr, 1)  // Implementation of the TPC pedestal and noise calibration
};

#endif

