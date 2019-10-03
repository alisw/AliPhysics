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
#include <TArrayD.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

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
  Bool_t GlobalQualityAssurance(AliAODEvent *aAODevent);
  Bool_t TrackSelection(AliAODTrack *aTrack); 
  virtual void CalculateQvectors();
  virtual void Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10);
  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 

  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetBoolMultCut(Bool_t top, Int_t otp){this->bMultCut = top; this->fSecondFilter=otp;} 

  void SetUpperLineCut(Float_t slope, Float_t axis)
  {this->fSlopeUpperLine = slope; this->fAxisUpperLine = axis; } 

   void SetLowerLineCut(Float_t slope, Float_t axis)
  {this->fSlopeLowerLine = slope; this->fAxisLowerLine = axis; } 
  
  void SetFilter(Int_t top){this->fMainFilter = top;} 

  void SetCentralityEstimator(Bool_t Esti){this->fCentralityfromVZero = Esti; }
 
  void SetVertexX(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnVertexX=Cut;  this->fMinVertexX=Min; this->fMaxVertexX=Max;}

  void SetVertexY(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnVertexY=Cut;  this->fMinVertexY=Min; this->fMaxVertexY=Max;}

  void SetVertexZ(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnVertexZ=Cut;  this->fMinVertexZ=Min; this->fMaxVertexZ=Max;}

   void SetEtaCut(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnEta=Cut;  this->fMinEtaCut=Min; this->fMaxEtaCut=Max;}

  void SetPtCut(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bCutOnPt=Cut;  this->fMinPtCut=Min; this->fMaxPtCut=Max;}

  void SetNumberTPCClusters(Bool_t Cut, Int_t Min)
  {this->bNumberTPCCluster=Cut;  this->fMinTPCCluster=Min; }

  void SetNumberITSClusters(Bool_t Cut, Int_t Min)
  {this->bNumberITSCluster=Cut;  this->fMinITSCluster=Min; }

  void SetChiSquareTPC(Bool_t Cut, Double_t Min, Double_t Max)
  {this->bChiSquareTPC=Cut;  this->fMinChiSquareTPC=Min; this->fMaxChiSquareTPC=Max;}

  void SetDCAz(Bool_t Cut, Double_t Max)
  {this->bDCAz=Cut;  this->fMaxDCAz=Max;}
 
  void SetMinNuPar(Int_t top){this->fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return this->fMinNumberPart;}

  void SetCorrSet1(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t i, Int_t j)
  {this->fNumber=Number; this->fh1=a; this->fh2=b; this->fh3=c; this->fh4=d; this->fh5=e; this->fh6=f; this->fh7=g; this->fh8=h; this->fh9=i; this->fh10=j;}

   void SetCorrSet2(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t i, Int_t j)
  {this->fNumberSecond=Number; this->fa1=a; this->fa2=b; this->fa3=c; this->fa4=d; this->fa5=e; this->fa6=f; this->fa7=g; this->fa8=h; this->fa9=i; this->fa10=j;}

  void SetRatioWeight(Bool_t top)
  {this->bUseRatioWeight=top;}

  void SetDenominatorMinValue(Double_t top) {this->fDenominatorMinValue=top; }

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
  TH1F *fPhiHistBeforeTrackSeletion;                // atrack->Phi() - Distribution before Track Selection
  TH1F *fEtaHistBeforeTrackSeletion;                // atrack->Eta() - Distribution before Track Selection
  TH1F *fTotalMultBeforeTrackSeletion;         // total number of Multiplicity for a centrality before Track Selection
  TH1F *fMultiHistoBeforeTrackSeletion;             // multiplicity distribution before Track Selection
  TH1F *fPhiHistAfterTrackSeletion;                // atrack->Phi() - Distribution before Track Selection
  TH1F *fEtaHistAfterTrackSeletion;                // atrack->Eta() - Distribution before Track Selection
  TH1F *fTotalMultAfterTrackSeletion;         // total number of Multiplicity for a centrality before Track Selection
  TH1F *fMultiHistoAfterTrackSeletion;             // multiplicity distribution before Track Selection
  TH1F *fMultiHistoBeforeMultCut;             // multiplicity distribution before high multiplicity outlier removel
  TH1F *fTPCClustersBeforeCut;		//Number of TPC clusters before cut
  TH1F *fTPCClustersAfterCut;		//Number of TPC clustes after cut
  TH1F *fITSClustersBeforeCut;		//Number of ITS clusters before cut
  TH1F *fITSClustersAfterCut;		//Number of ITS clusters after cut
  TH1F *fChiSquareTPCBeforeCut;		//Chi Square TPC before cut
  TH1F *fChiSquareTPCAfterCut;		//Chi Square TPC after cut
  TH1F *fDCAzBeforeCut;			//DCAz before cut
  TH1F *fDCAzAfterCut;			//DCAz after cut

  //2.) SelectionCuts
  Bool_t bMultCut;
  Int_t fMainFilter;           //for main filter selection (default: Hypbrid)
  Int_t fSecondFilter;           //for filter selection (default: global)
  Float_t fSlopeUpperLine;           //slope of the upper line for multiplicity cut
  Float_t fAxisUpperLine;           //axis intercept of the upper line for multiplicity cut
  Float_t fSlopeLowerLine;           //slope of the lower line for multiplicity cut
  Float_t fAxisLowerLine;           //axis intercept of the lower line for multiplicity cut
  
    //Global
  Float_t fMinCentrality;        // min centrality (default 0.)
  Float_t fMaxCentrality;        // max centrality (default 100.)
  Bool_t bCutOnVertexX;               // Bool to apply Vertex Cut in X (default kFALSE)
  Bool_t bCutOnVertexY;               // Bool to apply Vertex Cut in Y (default kFALSE)
  Bool_t bCutOnVertexZ;               // Bool to apply Vertex Cut in Z (default kFALSE)
  Double_t fMinVertexX;               // min vertex cut X (default -44)
  Double_t fMaxVertexX;               // max vertex cut X (default -44)
  Double_t fMinVertexY;               // min vertex cut Y (default -44)
  Double_t fMaxVertexY;               // max vertex cut Y (default -44)
  Double_t fMinVertexZ;               // min vertex cut Z (default -10 cm)
  Double_t fMaxVertexZ;               // max vertex cut Z (default +10 cm)
  TH1F *fVertexXBefore;               // Histogram Vertex X before vertex cut
  TH1F *fVertexXAfter;               // Histogram Vertex X after vertex cut
  TH1F *fVertexYBefore;               // Histogram Vertex Y before vertex cut
  TH1F *fVertexYAfter;               // Histogram Vertex Y after vertex cut
  TH1F *fVertexZBefore;               // Histogram Vertex Z before vertex cut
  TH1F *fVertexZAfter;               // Histogram Vertex Z after vertex cut
  Bool_t fCentralityfromVZero;	     // if kTRUE: Use V0 as centrality estimator, if kFALSE: SPD Cluster

    //Physics-Selection
  Bool_t bCutOnEta;               // Bool to apply eta cuts (default kTRUE)
  Bool_t bCutOnPt;               // Bool to apply pt cuts (default kTRUE)
  Bool_t bNumberTPCCluster;	//Bool to apply cuts on number of TPC clusters (default kTRUE)
  Bool_t bNumberITSCluster;	//Bool to apply cuts on number of ITS clusters (default kTRUE)
  Bool_t bChiSquareTPC;		//Bool to apply cuts on chi square TPC (default kTRUE)
  Bool_t bDCAz;			//Bool to apply cuts on DCAz (default kTRUE)
  Double_t fMinEtaCut;               // min eta cut (default -0.8)
  Double_t fMaxEtaCut;               // max eta cut (default 0.8)
  Double_t fMinPtCut;               // min pt cut (default 0.2)
  Double_t fMaxPtCut;               // max pt cut (default 5.0)
  Int_t fMinTPCCluster;		//Number of minimum TPC clusters (default 70)
  Int_t fMinITSCluster;		//Number of minimum ITS clusters (default 2)
  Double_t fMinChiSquareTPC;	//Minimal Chi Square TPC (default 0.1)
  Double_t fMaxChiSquareTPC;	//Maximal Chi Square TPC (default 4.0)
  Double_t fMaxDCAz;		//Maximal DCAz (default 2.4 cm)

  //3.) Variables for the correlation:
  Int_t fMaxCorrelator;          // maximum of correlation 
  TProfile *fRecursion[2][10];    //!  
  Bool_t bUseWeights;  
  
  Int_t fNumber;           //number of correlation first correlator
  Int_t fNumberSecond;           //number of correlation second correlator
  Int_t fMinNumberPart;           //minimal number of particles to do correlation
  Bool_t bUseRatioWeight;	//use number of combination weight for EbE Ratio (default kTRUE)
  Double_t fDenominatorMinValue;   //min value for the denominator in EbE Ratio (default 

  Int_t fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8, fh9, fh10;  //harmonics
  Int_t fa1, fa2, fa3, fa4, fa5, fa6, fa7, fa8, fa9, fa10;  //second set of harmonics


  Int_t fParticles;        // number of particles after all selections
  TArrayD *fAngles;              //! Azimuthal angles 
  TArrayD *fWeights;            //! Particle weights
  TArrayI *fBin;                   //! Bins for particle weight
  
  TComplex fQvector[61][11];       //! //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]

  // 4.) Final results:
   
  TProfile *fCentrality;         // final centrality result
  TProfile *fCentralitySecond;         // final centrality result for second harmonics 
  TProfile *fEvCentrality;         // final centrality result for event version
  TProfile *fCentralitySecondSquare; // final centrality result for second harmonics to the power of 2
  TProfile *fCentralitySecondSquareUnit; // final centrality result for second harmonics to the power of 2 with unit weights
  TProfile *fCov;         // Covariance term between first set of harmonics and second set of harmonics
  TProfile *fCovUnit;         // Covariance term between first set of harmonics and second set of harmonics
  TH1F *fCounterHistogram;       // for some checks
  TList *fFinalResultsList;      // list to hold all histograms with final results

  

  ClassDef(AliAnalysisTaskStudentsML,20);

};

//================================================================================================================

#endif
