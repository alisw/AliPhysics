/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************************** 
* template class for student projects         *
* author: Marcel Lesch (marcel.lesch@cern.ch) *
**********************************************/ 

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
#include <TArrayF.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"

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
 
  TComplex Q(Int_t n, Int_t p);
  TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip);
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();
	

  // 2.) Methods called in UserExec(Option_t *):
  // ...
  //add all except Cosmetics
  virtual void Cosmetics();
  virtual void CalculateQvectors();
  virtual void Correlation();
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

  /*void SetHolder(Int_t const maxcorrelators)
  {
   this->fMaxCorrelator = maxcorrelators; 
  };*/

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


  TH1F *fMultiHisto;             // multiplicity histogram atrack->nTracks

  //2.) Variables for the correlation:
  Int_t fMaxCorrelator;          // maximum of correlation 
  TProfile *fRecursion[2][8];    //! //how can i set the 8 as fMaxCorrelator?????? 
  Bool_t bUseWeights; 

  const Int_t kNumber;           //number of correlation

  const Int_t kh1, kh2, kh3, kh4, kh5, kh6, kh7, kh8;  //harmonics
   
  const Int_t kSum; 
  const Int_t kMaxHarmonic; 
  const Int_t kMaxPower; 
  Int_t fParticles;
  Float_t fCentral;
  Float_t fMinCentrality;        // min centrality
  Float_t fMaxCentrality;        // max centrality
  TArrayD *fAngles;              //! Azimuthal angles 
  TArrayD *fWeights;            //! Particle weights
  TArrayI *fBin;                   //! Bins for particle weight
  TF1 *func1;
  
  TComplex Qvector[17][9];       //! //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]

  // 3.) Final results:
   
  TProfile *fCentrality;         // final centrality result
  TH1F *fCounterHistogram;       // for some checks
  TList *fFinalResultsList;      // list to hold all histograms with final results

  

  ClassDef(AliAnalysisTaskStudentsML,4);

};

//================================================================================================================

#endif










