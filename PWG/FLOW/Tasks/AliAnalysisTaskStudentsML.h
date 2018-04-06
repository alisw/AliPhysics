/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* template class for student projects * 
**************************************/ 

#ifndef ALIANALYSISTASKSTUDENTSML_H
#define ALIANALYSISTASKSTUDENTSML_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TComplex.h"

//================================================================================================================

class AliAnalysisTaskStudentsML : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskStudentsML();
  AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskStudentsML(); 
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
  //add all except Cosmetics
  
  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 
  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbins = nbins;
   this->fMinBin = min;
   this->fMaxBin = max;
  };

 private:
  AliAnalysisTaskStudentsML(const AliAnalysisTaskStudentsML& aatmpf);
  AliAnalysisTaskStudentsML& operator=(const AliAnalysisTaskStudentsML& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TH1F *fPtHist;                 // atrack->Pt()
  Int_t fNbins;                  // number of bins
  Float_t fMinBin;               // min bin
  Float_t fMaxBin;               // min bin 
  TH1F *fPhiHist;                // atrack->Phi()
  TH1F *fEtaHist;                // atrack->Eta()
  TH1F *fMultiHist;		 // multiplicity histogram atrack->nTracks

  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  Int_t fMaxCorrelator;

  ClassDef(AliAnalysisTaskStudentsML,1);

};

//================================================================================================================

#endif











