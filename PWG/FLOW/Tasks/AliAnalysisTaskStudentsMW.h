/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/************************************** 
* template class for student projects * 
**************************************/ 

#ifndef ALIANALYSISTASKSTUDENTSMW_H
#define ALIANALYSISTASKSTUDENTSMW_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TComplex.h"
#include <TArrayD.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

//================================================================================================================

class AliAnalysisTaskStudentsMW : public AliAnalysisTaskSE{
 public:
  
  AliAnalysisTaskStudentsMW();
  AliAnalysisTaskStudentsMW(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskStudentsMW(); 
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
  virtual void Cosmetics();
  Bool_t TrackSelection(AliAODTrack *aTrack); 
  virtual void CalculateQvectors();
  virtual void Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8);

  // 3.) Setters and getters:
   void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 

  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}
  
  void SetFilter(Int_t top){this->fMainFilter = top;} 

  void SetVertexZ(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnVertexZ=Cut;  this->fMinVertexZ=Min; this->fMaxVertexZ=Max;}

   void SetEtaCut(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnEta=Cut;  this->fMinEtaCut=Min; this->fMaxEtaCut=Max;}

  void SetPtCut(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnPt=Cut;  this->fMinPtCut=Min; this->fMaxPtCut=Max;}
 
  void SetMinNuPar(Int_t top){this->fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return this->fMinNumberPart;}

  void SetCorrSet1(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h)
  {this->fNumber=Number; this->fh1=a; this->fh2=b; this->fh3=c; this->fh4=d; this->fh5=e; this->fh6=f; this->fh7=g; this->fh8=h;}

   void SetCorrSet2(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h)
  {this->fNumberSecond=Number; this->fa1=a; this->fa2=b; this->fa3=c; this->fa4=d; this->fa5=e; this->fa6=f; this->fa7=g; this->fa8=h;}

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

 private:
  AliAnalysisTaskStudentsMW(const AliAnalysisTaskStudentsMW& aatmpf);
  AliAnalysisTaskStudentsMW& operator=(const AliAnalysisTaskStudentsMW& aatmpf);
  
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

  TH1F *fPhiHistBeforeTrackSeletion;                // atrack->Phi() - Distribution before Track Selection
  TH1F *fEtaHistBeforeTrackSeletion;                // atrack->Eta() - Distribution before Track Selection
  TH1F *fTotalMultBeforeTrackSeletion;         // total number of Multiplicity for a centrality before Track Selection
  TH1F *fMultiHistoBeforeTrackSeletion;             // multiplicity distribution before Track Selection
  TH1F *fPhiHistAfterTrackSeletion;                // atrack->Phi() - Distribution before Track Selection
  TH1F *fEtaHistAfterTrackSeletion;                // atrack->Eta() - Distribution before Track Selection
  TH1F *fTotalMultAfterTrackSeletion;         // total number of Multiplicity for a centrality before Track Selection
  TH1F *fMultiHistoAfterTrackSeletion;             // multiplicity distribution before Track Selection
  TH1F *fMultiHistoBeforeMultCut;             // multiplicity distribution before high multiplicity 
  
  //2.) SelectionCuts
  Int_t fMainFilter;           //for main filter selection (default: Hypbrid)
  
    //Global
  Bool_t bCutOnVertexZ;               // Bool to apply Vertex Cut in Z (default kFALSE)
  Double_t fMinVertexZ;               // min vertex cut Z (default -10 cm)
  Double_t fMaxVertexZ;               // max vertex cut Z (default +10 cm)
  TH1F *fVertexZBefore;               // Histogram Vertex Z before vertex cut
  TH1F *fVertexZAfter;               // Histogram Vertex Z after vertex cut

    //Physics-Selection
  Bool_t bCutOnEta;               // Bool to apply eta cuts (default kTRUE)
  Bool_t bCutOnPt;               // Bool to apply pt cuts (default kTRUE)
  Double_t fMinEtaCut;               // min eta cut (default -0.8)
  Double_t fMaxEtaCut;               // max eta cut (default 0.8)
  Double_t fMinPtCut;               // min pt cut (default 0.2)
  Double_t fMaxPtCut;               // max pt cut (default 5.0)
 

  //3.) Variables for the correlation:
  Int_t fMaxCorrelator;          // maximum of correlation 
  TProfile *fRecursion[2][8];    //!  
  Bool_t bUseWeights; 

  
  Int_t fNumber;           //number of correlation first correlator
  Int_t fNumberSecond;           //number of correlation second correlator
  Int_t fMinNumberPart;           //minimal number of particles to do correlation

  Int_t fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8;  //harmonics
  Int_t fa1, fa2, fa3, fa4, fa5, fa6, fa7, fa8;  //second set of harmonics


  Int_t fParticles;        // number of particles after all selections
  TArrayD *fAngles;              //! Azimuthal angles 
  TArrayD *fWeights;            //! Particle weights
  TArrayI *fBin;                   //! Bins for particle weight
  
  TComplex fQvector[49][9];       //! //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]

  // 4.) Final results:
   
  TProfile *fCentralityres;         // final centrality result
  TProfile *fCentralitySecond;         // final centrality result for second harmonics 
  TProfile *fCentralitySecondSquare; // final centrality result for second harmonics to the power of 2
  TProfile *fCov;         // Covariance term between first set of harmonics and second set of harmonics
  TH1F *fCounterHistogram;       // for some checks
  TList *fFinalResultsList;      // list to hold all histograms with final results
  
   
  // new end

  ClassDef(AliAnalysisTaskStudentsMW,7);

};

//================================================================================================================

#endif

