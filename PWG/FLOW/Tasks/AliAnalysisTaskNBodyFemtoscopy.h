/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* Femtoscopy task for N bodies * 
Laura Serksnyte
laura.serksnyte@cern.ch
**************************************/ 

#ifndef ALIANALYSISTASKNBODYFEMTOSCOPY_H
#define ALIANALYSISTASKNBODYFEMTOSCOPY_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TString.h"

//================================================================================================================

class AliAnalysisTaskNBodyFemtoscopy : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskNBodyFemtoscopy();
  AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskNBodyFemtoscopy(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0.) Methods called in the constructor:
  virtual void InitializeArrays();
 
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *):
  // ...

  // 3.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetPtBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsPt = nbins;
   this->fMinBinPt = min;
   this->fMaxBinPt = max;
  };

  void SetCentralityEstimator(TString estimatorName)
  {
   this->fCentralityEstimator = estimatorName;
  };

  void SetCentralityBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbinsCentrality = nbins;
   this->fMinCentrality = min;
   this->fMaxCentrality = max;
  };


   private:
  AliAnalysisTaskNBodyFemtoscopy(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  AliAnalysisTaskNBodyFemtoscopy& operator=(const AliAnalysisTaskNBodyFemtoscopy& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TH1F *fPtHist;                 // atrack->Pt()
  TH1F *fCentralityHist;         // ams->GetMultiplicityPercentile()
  Int_t fNbinsPt;                // number of bins
  Float_t fMinBinPt;             // min bin
  Float_t fMaxBinPt;             // min bin
  TString fCentralityEstimator;  // centrality estimator choice
  Float_t fNbinsCentrality;      // cenrality bin number
  Float_t fMinCentrality;        // min centrality
  Float_t fMaxCentrality;        // max centrality
  
  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskNBodyFemtoscopy,2);

};

//================================================================================================================

#endif




