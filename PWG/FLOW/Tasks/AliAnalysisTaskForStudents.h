/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* template class for student projects * 
**************************************/ 

#ifndef ALIANALYSISTASKFORSTUDENTS_H
#define ALIANALYSISTASKFORSTUDENTS_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TGraphErrors.h"
#include "TComplex.h"

//================================================================================================================

class AliAnalysisTaskForStudents : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskForStudents();
  AliAnalysisTaskForStudents(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskForStudents(); 
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

  void SetNumValue(Int_t const dummy){this->num = dummy;};
  //TList* GetFinalResultsList() const {return this->fFinalResultsList;}


  void SetBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNbins = nbins;
   this->fMinBin = min;
   this->fMaxBin = max;
  };

  void SetCentralityBinning(Int_t const nbins, Float_t min, Float_t max)
  {
   this->fNCentralityBins = nbins;
   this->fMinCentrality = min;
   this->fMaxCentrality = max;
  };

  void Setylimit( Double_t ymin, Double_t ymax)
  {
   this->fymin = ymin;
   this->fymax = ymax;
  };

 private:
  AliAnalysisTaskForStudents(const AliAnalysisTaskForStudents& aatmpf);
  AliAnalysisTaskForStudents& operator=(const AliAnalysisTaskForStudents& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms:  
  TList *fControlHistogramsList; // list to hold all control histograms
  TH1F *fPtHist;                 // atrack->Pt()
  Int_t fNbins;                  // number of bins
  Float_t fMinBin;               // min bin
  Float_t fMaxBin;               // min bin
 
  TH1F *fCentralityHist;         // ams->GetMultiplicityPercentile("V0M")
  Int_t fNCentralityBins;        // number of centrality bins
  Float_t fMinCentrality;        // min centrality
  Float_t fMaxCentrality;        // max centrality

  TH1F *fMultHist;  
  TH1F *fPhiHist;  
  TH1F *fEtaHist;  

  TH1F *fhf;
  TH1F *fhf2;
  TLegend *leg_hist;

  

  // 2.) Final results:
  TList *fFinalResultsList; // list to hold all histograms with final results

  // new start
  Int_t  num, fhr2bin, fsize, fbinnum, fcounter1, fcounter2, fcounter3, fmult;
  Double_t c2, fhr2min, fhr2max, fper2, fper3, fper4, fvarb, fmu, fymin, fymax, fc4, fc3, fc22, fc24, fc23;
  
  Double_t fsc4[9], fcentral[9], fyerr[9]; // [0] = x, [1] = y

  TProfile *fhr2;
  TProfile *fhr3;
  TProfile *fhrc22;
  TProfile *fhrc23;
  TProfile *fhrc24;
  TProfile *fhrc3;
  TProfile *fhrc4;


  TGraphErrors *fgr;
  
  TComplex fq6, fq4, fq, fq2, fqq, fq5, fq3, fq1, fq3x, fqq3;
  
   
  // new end

  ClassDef(AliAnalysisTaskForStudents,3);

};

//================================================================================================================

#endif











