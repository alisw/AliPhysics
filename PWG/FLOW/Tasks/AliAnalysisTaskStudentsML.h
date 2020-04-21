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
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliHeader.h"
#include "TExMap.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TFile.h"
#include "TSystem.h"

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
  virtual void Cosmetics();	

  // 2.) Methods called in UserExec(Option_t *) or subsquent Methods:
  virtual void PhysicsAnalysis(AliAODEvent *aAODEvent);
  virtual void GetKineDist(AliAODEvent *aAODEve, AliMCEvent *aMCEve);

  //a) Methods used to assure global quality and track selection
  Bool_t GlobalQualityAssurance(AliAODEvent *aAODevent);
  Bool_t TrackSelection(AliAODTrack *aTrack); 
  Bool_t GlobalQualityAssurance(AliMCEvent *aMCKineEvent);
  Bool_t TrackSelection(AliAODMCParticle *aMCtrack);

  //b) Method used to obtain particle weights
  virtual void CalculateWeight(Int_t RunNumber, Double_t* weights, Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta);
  Int_t GetRunIndex(Int_t runNumber);
  //c) Methods used for Multi-Particle Correlation Techniques
  virtual void CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights);
  TComplex CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B);
  TComplex Q(Int_t n, Int_t p);
  TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip);
  virtual void Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Double_t* Correlation_Angle, Int_t Correlation_Mult, Double_t* Correlation_Weight);

  virtual void MainTask(Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array);
  virtual void MixedParticle(Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B);

  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters:
  void SetControlHistogramsList(TList* const chl) {this->fControlHistogramsList = chl;};
  TList* GetControlHistogramsList() const {return this->fControlHistogramsList;} 

  void SetFinalResultsList(TList* const frl) {this->fFinalResultsList = frl;};
  TList* GetFinalResultsList() const {return this->fFinalResultsList;}

  void SetDoAnalysis(Bool_t DoAnalysis, Bool_t RicoKineTable){this->bDoAnalysis = DoAnalysis; this->bUseRecoKineTable = RicoKineTable;}

  void SetBoolMultCut(Bool_t top, Int_t otp){this->bMultCut = top; this->fSecondFilter=otp;} 

  void SetUpperLineCut(Float_t slope, Float_t axis)
  {this->fSlopeUpperLine = slope; this->fAxisUpperLine = axis; } 

   void SetLowerLineCut(Float_t slope, Float_t axis)
  {this->fSlopeLowerLine = slope; this->fAxisLowerLine = axis; } 
  
  void SetFilter(Int_t top){this->fMainFilter = top;} 

  void SetCentralityEstimator(Bool_t Esti){this->fCentralityfromVZero = Esti; }
 
  void SetUseWeights(Bool_t Weights, Bool_t Phi, Bool_t Pt, Bool_t Eta)
  {this->bUseWeights = Weights; this->bUsePhiWeights = Phi; this->bUsePtWeights = Pt; this->bUseEtaWeights = Eta;}

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
 
  void SetDCAxy(Bool_t Cut, Double_t Max)
  {this->bDCAxy=Cut;  this->fMaxDCAxy=Max;}

  void SetMinNuPar(Int_t top){this->fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return this->fMinNumberPart;}

  void SetCorrSet1(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t i, Int_t j, Int_t k, Int_t l)
  {this->fNumber=Number; this->fh1=a; this->fh2=b; this->fh3=c; this->fh4=d; this->fh5=e; this->fh6=f; this->fh7=g; this->fh8=h; this->fh9=i; this->fh10=j; this->fh11=k; this->fh12=l;}

   void SetCorrSet2(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t i, Int_t j, Int_t k, Int_t l)
  {this->fNumberSecond=Number; this->fa1=a; this->fa2=b; this->fa3=c; this->fa4=d; this->fa5=e; this->fa6=f; this->fa7=g; this->fa8=h; this->fa9=i; this->fa10=j;this->fa11=k; this->fa12=l;}

void SetCorrSet3(Bool_t booly, Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g, Int_t h, Int_t i, Int_t j, Int_t k, Int_t l)
  {this->bDoThirdCorrelation=booly; this->fNumberThird=Number; this->fb1=a; this->fb2=b; this->fb3=c; this->fb4=d; this->fb5=e; this->fb6=f; this->fb7=g; this->fb8=h; this->fb9=i; this->fb10=j;this->fb11=k; this->fb12=l;}

  void SetDoEbE(Bool_t top){this->bDoEbERatio=top;}

  void SetRatioWeight(Bool_t top){this->bUseRatioWeight=top;}

  void SetDenominatorMinValue(Double_t top) {this->fDenominatorMinValue=top; }

  void SetMixed(Bool_t top, Int_t nop, Bool_t DifferentCharge, Bool_t PositiveCharge)
  {this->bDoMixed = top; this->fMixedHarmonic = nop; this->bDifferentCharge = DifferentCharge; 
	this->bSetSameChargePositiv = PositiveCharge;	}

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

  void SetListOfRuns(TString dataPeriod);

  void SetInputParticleWeights(TString fileWeight);

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
  TList *fControlHistogramsList; 			// list to hold all control histograms
  Int_t fNbins;                  			// number of bins
  Float_t fMinBin;               			// min bin
  Float_t fMaxBin;              			// min bin 
 
  TH1F *fPTHistogram[4]; 		//! 0: P_t Before Track Selection, 1: P_t After Track Selection, 2: P_t Before Track Selection (second), 3: P_t After Track Selection (second);
  TH1F *fPhiHistogram[4]; 		//! 0: Phi Before Track Selection, 1: Phi After Track Selection, 2: Phi Before Track Selection (second), 3: Phi After Track Selection (second);
  TH1F *fEtaHistogram[4]; 		//! 0: Eta Before Track Selection, 1: Eta After Track Selection, 2: Eta Before Track Selection (second), 3: Eta After Track Selection (second);
  TH1F *fMultHistogram[4]; 		//! 0: Multiplicity before HMO removel 1: Multiplicity Before Track Selection, 2: Mult. After Track Selection 3: Mult. After Track Selection (second);
  TH1F *fTPCClustersHistogram[2]; 	//! 0: TPC Clusters Before Track Selection, 1: TPC Clusters After Track Selection
  TH1F *fITSClustersHistogram[2]; 	//! 0: ITS Clusters Before Track Selection, 1: ITS Clusters After Track Selection
  TH1F *fChiSquareTPCHistogram[2]; 	//! 0: ChiSquare TPC Before Track Selection, 1: ChiSquare TPC After Track Selection
  TH1F *fDCAzHistogram[2]; 		//! 0: DCAz Before Track Selection, 1: DCAz After Track Selection
  TH1F *fDCAxyHistogram[2]; 		//! 0: DCAxy Before Track Selection, 1: DCAxy After Track Selection
  TH1F *fCentralityHistogram[2]; 	//! 0: Centrality Before Corresponding, 1: Centrality After Corresponding Cut

  //2.) SelectionCuts
  Bool_t bDoAnalysis;			// if kTRUE: Run of AODs (real Data or MC on Recon level) Does Analysis, if kFALSE: Get Weights
  Bool_t bUseRecoKineTable;		// (Necessary if bDoAnalysis = kFALSE) if kTRUE: use the kine-reco mapping table (for kine level)
  Bool_t bMultCut;			// Cut to remove HMO's
  Int_t fMainFilter;           		// for main filter selection (default: Hypbrid)
  Int_t fSecondFilter;           	// for filter selection (default: global)
  Float_t fSlopeUpperLine;          	// slope of the upper line for multiplicity cut
  Float_t fAxisUpperLine;          	// axis intercept of the upper line for multiplicity cut
  Float_t fSlopeLowerLine;          	// slope of the lower line for multiplicity cut
  Float_t fAxisLowerLine;           	// axis intercept of the lower line for multiplicity cut
  
    //Global
  Float_t fMinCentrality;        	// min centrality (default 0.)
  Float_t fMaxCentrality;        	// max centrality (default 100.)
  Bool_t bCutOnVertexX;               	// Bool to apply Vertex Cut in X (default kFALSE)
  Bool_t bCutOnVertexY;               	// Bool to apply Vertex Cut in Y (default kFALSE)
  Bool_t bCutOnVertexZ;               	// Bool to apply Vertex Cut in Z (default kFALSE)
  Double_t fMinVertexX;               	// min vertex cut X (default -44)
  Double_t fMaxVertexX;               	// max vertex cut X (default -44)
  Double_t fMinVertexY;               	// min vertex cut Y (default -44)
  Double_t fMaxVertexY;               	// max vertex cut Y (default -44)
  Double_t fMinVertexZ;               	// min vertex cut Z (default -10 cm)
  Double_t fMaxVertexZ;               	// max vertex cut Z (default +10 cm)
  TH1F *fVertexXHistogram[2];               	//! 0: Vertex X Before Corresponding, 1: Vertex X After Corresponding Cut
  TH1F *fVertexYHistogram[2];               	//! 0: Vertex Y Before Corresponding, 1: Vertex Y After Corresponding Cut
  TH1F *fVertexZHistogram[2];               	//! 0: Vertex Z Before Corresponding, 1: Vertex Z After Corresponding Cut

  Bool_t fCentralityfromVZero;	     	// if kTRUE: Use V0 as centrality estimator, if kFALSE: SPD Cluster

    //Physics-Selection
  Bool_t bCutOnEta;               	// Bool to apply eta cuts (default kTRUE)
  Bool_t bCutOnPt;               	// Bool to apply pt cuts (default kTRUE)
  Bool_t bNumberTPCCluster;		// Bool to apply cuts on number of TPC clusters (default kTRUE)
  Bool_t bNumberITSCluster;		// Bool to apply cuts on number of ITS clusters (default kTRUE)
  Bool_t bChiSquareTPC;			// Bool to apply cuts on chi square TPC (default kTRUE)
  Bool_t bDCAz;				// Bool to apply cuts on DCAz (default kTRUE)
  Bool_t bDCAxy;			// Bool to apply cuts on DCAxy (default kTRUE)
  Double_t fMinEtaCut;               	// min eta cut (default -0.8)
  Double_t fMaxEtaCut;               	// max eta cut (default 0.8)
  Double_t fMinPtCut;               	// min pt cut (default 0.2)
  Double_t fMaxPtCut;               	// max pt cut (default 5.0)
  Int_t fMinTPCCluster;			// Number of minimum TPC clusters (default 70)
  Int_t fMinITSCluster;			// Number of minimum ITS clusters (default 2)
  Double_t fMinChiSquareTPC;		// Minimal Chi Square TPC (default 0.1)
  Double_t fMaxChiSquareTPC;		// Maximal Chi Square TPC (default 4.0)
  Double_t fMaxDCAz;			// Maximal DCAz (default 3.2 cm)
  Double_t fMaxDCAxy;			// Maximal DCAxy (default 2.4 cm)

  //Weights
  Bool_t bUseWeights; 
  Bool_t bUsePtWeights;
  Bool_t bUsePhiWeights; 
  Bool_t bUseEtaWeights;
  Int_t fNumberRuns;            // Number of runs in the dataset.
  Int_t fListRuns[90];          // List of runs in the dataset.
  TH1F *fHistoPtWeight[90];     // Histograms with the pT-weights for each run.
  TH1F *fHistoEtaWeight[90];    // Histograms with the eta-weights for each run.
  TH1F *fHistoPhiWeight[90];    // Histograms with the phi-weights for each run.

  //3.) Variables for the correlation:
  Int_t fMaxCorrelator;          	// maximum of correlation 
  TProfile *fRecursion[2][12];    	//!   
  
  Int_t fNumber;           		// Number of correlation first correlator
  Int_t fNumberSecond;          	// Number of correlation second correlator
  Int_t fNumberThird;           	// Number of correlation third correlator
  Bool_t bDoThirdCorrelation;   	// if kTRUE: do the third correlation (default kFALSE)
  Int_t fMinNumberPart;           	// Minimal number of particles to do correlation
  Bool_t bUseRatioWeight;		// use number of combination weight for EbE Ratio (default kTRUE)
  Double_t fDenominatorMinValue;   	// min value for the denominator in EbE Ratio (default 

  Int_t fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8, fh9, fh10, fh11, fh12;  //harmonics
  Int_t fa1, fa2, fa3, fa4, fa5, fa6, fa7, fa8, fa9, fa10, fa11, fa12;  //second set of harmonics
  Int_t fb1, fb2, fb3, fb4, fb5, fb6, fb7, fb8, fb9, fb10, fb11, fb12;  //third set of harmonics
  
  TComplex fQvector[97][13];       	//! //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]

  // 4.) Final results:
   
  TProfile *fCentrality;         	// final centrality result
  TProfile *fCentralitySecond;         	// final centrality result for second harmonics 
  TProfile *fCentralityThird;         	// final centrality result for third harmonics 
  TProfile *fEvCentrality;         	// final centrality result for event version
  Bool_t bDoEbERatio;		 	// if kTRUE: Do the EbE ratio, Default: kFALSE
  TProfile *fMixedParticleHarmonics; 	// Stores output for special mixed particle analysis
  Bool_t bDoMixed;		 	// if kTRUE: Do special mixed particle analysis, default kFALSE (MainTask)
  Bool_t bDifferentCharge; 	 	// used in DoMixed: if kTRUE mixed particle analysis between positiv and negativ
				 	//		    if kFALSE mixed particle analysis between same charge 
				 	//		    (only positiv or only negativ particles)
				 	// Default kTRUE
  Bool_t bSetSameChargePositiv;   	// used if bDifferentCharge: if kTRUE use positiv, if kFALSE use negative (default kTRUE)
  Int_t fMixedHarmonic;			// Harmonic of special mixed particle analysis
  TH1F *fCounterHistogram;       	// for some checks
  TList *fFinalResultsList;      	// list to hold all histograms with final results

  

  ClassDef(AliAnalysisTaskStudentsML,28);

};

//================================================================================================================

#endif
