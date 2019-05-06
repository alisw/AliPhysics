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

  void SetBoolMultCut(Bool_t top){this->bBruteMultCut = top;} 
  Bool_t GetBoolMultCut() const {return this->bBruteMultCut;}

  void SetMultCut(Int_t top){this->fMultCut = top;} 
  Int_t GetMultCut() const {return this->fMultCut;}

  void SetFilter(Int_t top){this->fFilter = top;} 
  Int_t GetFilter() const {return this->fFilter;}

  void SetCorrNumber(Int_t top){this->fNumber = top;} 
  Int_t GetCorrNumber() const {return this->fNumber;}

  void SetMinNuPar(Int_t top){this->fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return this->fMinNumberPart;}

  void SetCorrSet1(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h)
  {this->fh1=a; this->fh2=b; this->fh3=c; this->fh4=d; this->fh5=e; this->fh6=f; this->fh7=g; this->fh8=h;}

   void SetCorrSet2(Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h)
  {this->fa1=a; this->fa2=b; this->fa3=c; this->fa4=d; this->fa5=e; this->fa6=f; this->fa7=g; this->fa8=h;}

  void SetMinCent(Float_t top){this->fMinCentrality = top;} 
  Float_t GetMinCent() const {return this->fMinCentrality;}

  void SetMaxCent(Float_t top){this->fMaxCentrality = top;} 
  Float_t GetMaxCent() const {return this->fMaxCentrality;}
  

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
  TH1F *fMultPreCut;         // Multiplicity before brute cut
  TH1F *fMultPostCut;         // Multiplicity after brute cut
  TH1F *fMultiHisto;             // multiplicity histogram atrack->nTracks

  //2.) Variables for the correlation:
  Int_t fMaxCorrelator;          // maximum of correlation 
  TProfile *fRecursion[2][8];    //!  
  TProfile *fRecursionSecond[2][8];    //!
  Bool_t bUseWeights; 

  Bool_t bBruteMultCut;
  Int_t fMultCut;
  Int_t fFilter;           //for filter selection
  Int_t fNumber;           //number of correlation
  Int_t fMinNumberPart;           //minimal number of particles to do correlation

  Int_t fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8;  //harmonics
  Int_t fa1, fa2, fa3, fa4, fa5, fa6, fa7, fa8;  //second set of harmonics

   
  const Int_t kSum; 
  const Int_t kMaxHarmonic; 
  const Int_t kMaxPower; 
  Int_t fParticles;
  Float_t fMinCentrality;        // min centrality
  Float_t fMaxCentrality;        // max centrality
  TArrayD *fAngles;              //! Azimuthal angles 
  TArrayD *fWeights;            //! Particle weights
  TArrayI *fBin;                   //! Bins for particle weight
  
  TComplex Qvector[17][9];       //! //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]

  // 3.) Final results:
   
  TProfile *fCentrality;         // final centrality result
  TProfile *fCentralitySecond;         // final centrality result for second harmonics 
  TProfile *fEvCentrality;         // final centrality result for event version
  TH1F *fCounterHistogram;       // for some checks
  TList *fFinalResultsList;      // list to hold all histograms with final results

  

  ClassDef(AliAnalysisTaskStudentsML,9);

};

//================================================================================================================

#endif










