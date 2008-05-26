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
  Bool_t Process();
  Bool_t ReadEntry(Int_t evt);
  Bool_t Accept(AliTPCclusterKr *cluster);
  Bool_t Update(AliTPCclusterKr *cluster);
  TH3F*  CreateHisto(Int_t chamber);

  const TObjArray* GetHistoKrArray () {return &fHistoKrArray;}  // get calibration object
  TH3F* GetHistoKr(Int_t sector) const;                         // get refernce histogram

  Bool_t IsCSide(Int_t chamber);
  Bool_t IsIROC(Int_t chamber);

  Bool_t Terminate();
  void Analyse();
  static TH1D* ProjectHisto(TH3F* histo3D, const char* name = "_pz", Int_t firstxbin = 0, Int_t lastxbin = 0, Int_t firstybin = 0, Int_t lastybin = 0);

  void SetInputChain(TChain *inChain) {fTree = (TTree*) inChain;}
  void SetOutputHisto(Bool_t bHisto = kTRUE) {bOutputHisto = bHisto;}

  void SetASide(Bool_t bA = kTRUE) {bASide = bA;} // fill histo only A TPC side
  void SetBSide(Bool_t bC = kTRUE) {bCSide = bC;} // fill histo only C TPC side

private:

  Bool_t bOutputHisto;        //!
  Bool_t bASide;              //!
  Bool_t bCSide;              //!

  TClonesArray *fClusters;    //!
  AliTPCclusterKr *fClustKr;  //!
  TTree *fTree;               //!
  TObjArray fHistoKrArray;    //  Calibration histograms for Kr distribution

public:
  ClassDef(AliTPCCalibKr, 1)  // Implementation of the TPC pedestal and noise calibration
};

#endif

