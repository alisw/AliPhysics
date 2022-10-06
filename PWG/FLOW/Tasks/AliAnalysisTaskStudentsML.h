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
  Bool_t TrackSelection(AliAODVertex *aPrimaryVertex, AliAODTrack *aTrack); 
  Bool_t GlobalQualityAssurance(Int_t CentBin, AliMCEvent *aMCKineEvent);
  Bool_t TrackSelection(AliAODMCParticle *aMCtrack);
  Int_t SelectCentrality(AliAODEvent *aAODevent);
  virtual void FisherYatesRandomizing(Int_t Mult, Int_t *RandomIndex); 

  //b) Method used to obtain particle weights
  virtual void CalculateWeight(Int_t CentBin, Int_t RunNumber, Double_t* weights, Int_t Multi, Double_t* angles, Double_t* pt, Double_t* eta);
  Int_t GetRunIndex(Int_t runNumber);
  //c) Methods used for Multi-Particle Correlation Techniques
  virtual void CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights);
  TComplex CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B);
  TComplex Q(Int_t n, Int_t p);
  TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip);
  virtual void Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Int_t h13, Int_t h14, Double_t *Correlation_Data);

  virtual void MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array);
  virtual void MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B);

  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters:
  void SetDoAnalysis(Bool_t DoAnalysis){this->bDoAnalysis = DoAnalysis; }

  void SetKineSpecifics(Bool_t RicoKineTable, Bool_t UseWeakSecondaries){this->bUseRecoKineTable = RicoKineTable; this->bUseWeakSecondaries = UseWeakSecondaries; } 

  void SetSaveAllQA(Bool_t SaveQA){this->bSaveAllQA=SaveQA;}

  void SetBoolMultCut(Bool_t top, Int_t otp){this->bMultCut = top; this->fSecondFilter=otp;} 

  void SetUpperLineCut(Float_t slope, Float_t axis)
  {this->fSlopeUpperLine = slope; this->fAxisUpperLine = axis; } 

   void SetLowerLineCut(Float_t slope, Float_t axis)
  {this->fSlopeLowerLine = slope; this->fAxisLowerLine = axis; } 
  
  void SetFilter(Int_t top){this->fMainFilter = top;} 

  void SetCentralityEstimator(TString Esti){this->fCentralityEstimator = Esti; }
 
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

  void SetChiSquareTPC(Bool_t Cut, Int_t ChooseMethod, Double_t Min, Double_t Max) 
  {this->bChiSquareTPC=Cut; this->fChooseChiSquareMethod=ChooseMethod; this->fMinChiSquareTPC=Min; this->fMaxChiSquareTPC=Max;}

  void SetDCAz(Bool_t Cut, Double_t Max)
  {this->bDCAz=Cut;  this->fMaxDCAz=Max;}
 
  void SetDCAxy(Bool_t Cut, Double_t Max)
  {this->bDCAxy=Cut;  this->fMaxDCAxy=Max;}

  void SetChargeCut(Bool_t Cut, Bool_t Pos)
  { this->bChargeCut=Cut; this->bChargePos=Pos; } 

  void SetFisherYates(Bool_t DoFY, Float_t CutOff)
  { this->bDoFisherYates=DoFY; this->fFisherYatesCutOff=CutOff; } 

  void SetMinNuPar(Int_t top){this->fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return this->fMinNumberPart;}

   void SetCorrSet1(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {this->fNumber=Number; this->fa1=a; this->fa2=b; this->fa3=c; this->fa4=d; this->fa5=e; this->fa6=f; this->fa7=g;}

  void SetCorrSet2( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { this->fNumberSecond=Number; this->fb1=a; this->fb2=b; this->fb3=c; this->fb4=d; this->fb5=e; this->fb6=f; this->fb7=g;}

  void SetCorrSet3( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { this->fNumberThird=Number; this->fd1=a; this->fd2=b; this->fd3=c; this->fd4=d; this->fd5=e; this->fd6=f; this->fd7=g;}

  void SetCorrSet4( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { this->fNumberFourth=Number; this->fe1=a; this->fe2=b; this->fe3=c; this->fe4=d; this->fe5=e; this->fe6=f; this->fe7=g;}

  void SetCorrSet5(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {this->fNumberFifth=Number; this->ff1=a; this->ff2=b; this->ff3=c; this->ff4=d; this->ff5=e; this->ff6=f; this->ff7=g;}

  void SetCorrSet6(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {this->fNumberSixth=Number; this->fg1=a; this->fg2=b; this->fg3=c; this->fg4=d; this->fg5=e; this->fg6=f; this->fg7=g;}

  void SetCorrSet7(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {this->fNumberSeventh=Number; this->fh1=a; this->fh2=b; this->fh3=c; this->fh4=d; this->fh5=e; this->fh6=f; this->fh7=g;}

  void SetCorrSet8(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {this->fNumberEighth=Number; this->fi1=a; this->fi2=b; this->fi3=c; this->fi4=d; this->fi5=e; this->fi6=f; this->fi7=g;}

  void SetMixed(Bool_t top, Int_t nop, Bool_t DifferentCharge, Bool_t PositiveCharge)
  {this->bDoMixed = top; this->fMixedHarmonic = nop; this->bDifferentCharge = DifferentCharge; 
	this->bSetSameChargePositiv = PositiveCharge;	}

  void SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4, Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10, Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16 )
 {this->fcent_0 = cen0; this->fcent_1 = cen1; this->fcent_2 = cen2; this->fcent_3 = cen3; this->fcent_4 = cen4; this->fcent_5 = cen5; this->fcent_6 = cen6; this->fcent_7 = cen7; this->fcent_8 = cen8; this->fcent_9 = cen9; this->fcent_10 = cen10; this->fcent_11 = cen11; this->fcent_12 = cen12; this->fcent_13 = cen13; this->fcent_14 = cen14; this->fcent_15 = cen15; this->fcent_16 = cen16;} 
 
  void SetListOfRuns(TString dataPeriod);

  void SetInputParticleWeights(TString fileWeight); //Set weight files that should be accessed

  void SetInitializeCentralityArray(); //Set Centrality array 


 private:
  AliAnalysisTaskStudentsML(const AliAnalysisTaskStudentsML& aatmpf);
  AliAnalysisTaskStudentsML& operator=(const AliAnalysisTaskStudentsML& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  // 1.) Control histograms: 
  TList *fCentralityList[16];		//! Will be one list per certain centrality bin. Up to 16 centraliy bins possible
  TList *fControlHistogramsList[16]; 	//! List to hold all control histograms for a specific centrality bin. Up to 16 centraliy bins possible
 
  TH1F *fPTHistogram[16][5]; 		//! 0: P_t Before Track Selection, 1: P_t After Track Selection, 2: P_t Before Track Selection (second), 3: P_t After Track Selection (second), 4: P_t After Track Selection Weighted; 
  TH1F *fPhiHistogram[16][5]; 		//! 0: Phi Before Track Selection, 1: Phi After Track Selection, 2: Phi Before Track Selection (second), 3: Phi After Track Selection (second), 4: Phi After Track Selection Weighted; 
  TH1F *fEtaHistogram[16][5]; 		//! 0: Eta Before Track Selection, 1: Eta After Track Selection, 2: Eta Before Track Selection (second), 3: Eta After Track Selection (second), 4: Eta After Track Selection Weighted; 
  TH1F *fMultHistogram[16][4]; 		//! 0: Multiplicity before HMO removel 1: Multiplicity Before Track Selection, 2: Mult. After Track Selection 3: Mult. After Track Selection (second);
  TH1F *fTPCClustersHistogram[16][2]; 	//! 0: TPC Clusters Before Track Selection, 1: TPC Clusters After Track Selection
  TH1F *fITSClustersHistogram[16][2]; 	//! 0: ITS Clusters Before Track Selection, 1: ITS Clusters After Track Selection
  TH1F *fChiSquareTPCHistogram[16][2]; 	//! 0: ChiSquare TPC Before Track Selection, 1: ChiSquare TPC After Track Selection
  TH1F *fDCAzHistogram[16][2]; 		//! 0: DCAz Before Track Selection, 1: DCAz After Track Selection
  TH1F *fDCAxyHistogram[16][2]; 	//! 0: DCAxy Before Track Selection, 1: DCAxy After Track Selection
  TH1I *fChargeHistogram[16][2]; 	//! 0: Charge Before Track Selection, 1: Charge After Track Selection 
  TH1F *fCentralityHistogram[16]; 	//! Centrality After Corresponding Cut
  TH1F *fCentralityHistogramBefore;     //! Centrality Histogram before Centrality selection

  //2.) SelectionCuts
  Bool_t bDoAnalysis;			// if kTRUE: Run of AODs (real Data or MC on Recon level) Does Analysis, if kFALSE: Get Weights
  Bool_t bUseRecoKineTable;		// (Necessary if bDoAnalysis = kFALSE) if kTRUE: use the kine-reco mapping table (for kine level)
  Bool_t bUseWeakSecondaries; 		// if kTRUE: use weak secondaries for obtaining weights, if kFALSE: only primaries
  Bool_t bSaveAllQA;			// if kTRUE: All Standard QA Histograms are saved (default kTRUE)
  Bool_t bMultCut;			// Cut to remove HMO's
  Int_t fMainFilter;           		// for main filter selection (default: Hypbrid)
  Int_t fSecondFilter;           	// for filter selection (default: global)
  Float_t fSlopeUpperLine;          	// slope of the upper line for multiplicity cut
  Float_t fAxisUpperLine;          	// axis intercept of the upper line for multiplicity cut
  Float_t fSlopeLowerLine;          	// slope of the lower line for multiplicity cut
  Float_t fAxisLowerLine;           	// axis intercept of the lower line for multiplicity cut

    //Centrality
  Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
  					//fcent_i holds the edge of a centrality bin
  Int_t fCentralityBins;		//will be set to 16, for at maximum 16 bins in case of centrality 0 to 80 in 5% steps. Less bins and different steps may be used
  Float_t fcentralityArray[17];		//will hold our edges of the centrality bins

    //Global
  Bool_t bCutOnVertexX;               	// Bool to apply Vertex Cut in X (default kFALSE)
  Bool_t bCutOnVertexY;               	// Bool to apply Vertex Cut in Y (default kFALSE)
  Bool_t bCutOnVertexZ;               	// Bool to apply Vertex Cut in Z (default kFALSE)
  Double_t fMinVertexX;               	// min vertex cut X (default -44)
  Double_t fMaxVertexX;               	// max vertex cut X (default -44)
  Double_t fMinVertexY;               	// min vertex cut Y (default -44)
  Double_t fMaxVertexY;               	// max vertex cut Y (default -44)
  Double_t fMinVertexZ;               	// min vertex cut Z (default -10 cm)
  Double_t fMaxVertexZ;               	// max vertex cut Z (default +10 cm)
  TH1F *fVertexXHistogram[16][2];       //! 0: Vertex X Before Corresponding, 1: Vertex X After Corresponding Cut
  TH1F *fVertexYHistogram[16][2];       //! 0: Vertex Y Before Corresponding, 1: Vertex Y After Corresponding Cut
  TH1F *fVertexZHistogram[16][2];       //! 0: Vertex Z Before Corresponding, 1: Vertex Z After Corresponding Cut

  TString fCentralityEstimator;	     	// Choose between: "V0M" for V0 as centrality estimator, or "CL1" for SPC Clusters
  //Physics-Selection
  Bool_t bCutOnEta;               	// Bool to apply eta cuts (default kTRUE)
  Bool_t bCutOnPt;               	// Bool to apply pt cuts (default kTRUE)
  Bool_t bNumberTPCCluster;		// Bool to apply cuts on number of TPC clusters (default kTRUE)
  Bool_t bNumberITSCluster;		// Bool to apply cuts on number of ITS clusters (default kTRUE)
  Bool_t bChiSquareTPC;			// Bool to apply cuts on chi square TPC (default kTRUE)
  Bool_t bDCAz;				// Bool to apply cuts on DCAz (default kTRUE)
  Bool_t bDCAxy;			// Bool to apply cuts on DCAxy (default kTRUE)
  Bool_t bChargeCut;			// Bool to apply cut on charge (default kFALSE) 
  Bool_t bChargePos;			// Bool to select only positive charge (if kTRUE) and negative if (kFALSE). Only relevant if bChargeCut==kTRUE 
  Double_t fMinEtaCut;               	// min eta cut (default -0.8)
  Double_t fMaxEtaCut;               	// max eta cut (default 0.8)
  Double_t fMinPtCut;               	// min pt cut (default 0.2)
  Double_t fMaxPtCut;               	// max pt cut (default 5.0)
  Int_t fMinTPCCluster;			// Number of minimum TPC clusters (default 70)
  Int_t fMinITSCluster;			// Number of minimum ITS clusters (default 2)
  Int_t fChooseChiSquareMethod;		// Choose how Chi Square is defined
  Double_t fMinChiSquareTPC;		// Minimal Chi Square TPC (default 0.1)
  Double_t fMaxChiSquareTPC;		// Maximal Chi Square TPC (default 4.0)
  Double_t fMaxDCAz;			// Maximal DCAz (default 3.2 cm)
  Double_t fMaxDCAxy;			// Maximal DCAxy (default 2.4 cm)

  //Fisher-Yales
  Bool_t bDoFisherYates;		//if kTRUE: Do Fisher Yates Mixing of phi, pt and eta arrays after track selection (default: kFALSE)
  Float_t fFisherYatesCutOff;		//How much percentage of the orginal particles are kept, e.g. if 0.7 only 70% of the current particles are kept for analysis

  //Weights
  Bool_t bUseWeights; 
  Bool_t bUsePtWeights;
  Bool_t bUsePhiWeights; 
  Bool_t bUseEtaWeights;
  Int_t fNumberRuns;           	    // Number of runs in the dataset.
  Int_t fListRuns[90];              // List of runs in the dataset.
  TH1F *fHistoPtWeight[16][90];     // Histograms with the pT-weights for each run.
  TH1F *fHistoEtaWeight[16][90];    // Histograms with the eta-weights for each run.
  TH1F *fHistoPhiWeight[16][90];    // Histograms with the phi-weights for each run.

  //3.) Variables for the correlation:
  Int_t fMaxCorrelator;          	// maximum of correlation   
  
  Int_t fNumber;           		// Number of correlation first correlator
  Int_t fNumberSecond;          	// Number of correlation second correlator
  Int_t fNumberThird;           	// Number of correlation third correlator
  Int_t fNumberFourth;			// Number of correlation fourth correlator
  Int_t fNumberFifth;			// Number of correlation fifth correlator
  Int_t fNumberSixth;			// Number of correlation sixth correlator
  Int_t fNumberSeventh;			// Number of correlation seventh correlator
  Int_t fNumberEighth;			// Number of correlation eigth correlator
  Int_t fMinNumberPart;           	// Minimal number of particles to do correlation

  Int_t fa1, fa2, fa3, fa4, fa5, fa6, fa7; //first set of harmonics
  Int_t fb1, fb2, fb3, fb4, fb5, fb6, fb7; //second set of harmonics
  Int_t fd1, fd2, fd3, fd4, fd5, fd6, fd7; //third set of harmonics
  Int_t fe1, fe2, fe3, fe4, fe5, fe6, fe7; //fourth set of harmonics
  Int_t ff1, ff2, ff3, ff4, ff5, ff6, ff7;  //sixth set of harmonics 
  Int_t fg1, fg2, fg3, fg4, fg5, fg6, fg7;  //seventh set of harmonics 
  Int_t fh1, fh2, fh3, fh4, fh5, fh6, fh7;  //harmonics
  Int_t fi1, fi2, fi3, fi4, fi5, fi6, fi7;  //eigth set of harmonics 
  
  TComplex fQvector[113][15];       	// //[fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1] 

  // 4.) Final results:
   
  TProfile *fResults[16];         	//! final centrality result
  TProfile *fCovResults[16];         	//! TProfile to store terms needed for Covariance 
  TProfile *fJoinedCovResults[16];      //! TProfile to store joined Covariance term calculated as one correlator <z> instead of product of two correlators <x*y> 
  TProfile *fMixedParticleHarmonics[16];//! Stores output for special mixed particle analysis
  Bool_t bDoMixed;		 	// if kTRUE: Do special mixed particle analysis, default kFALSE (MainTask)
  Bool_t bDifferentCharge; 	 	// used in DoMixed: if kTRUE mixed particle analysis between positiv and negativ
				 	//		    if kFALSE mixed particle analysis between same charge 
				 	//		    (only positiv or only negativ particles)
				 	// Default kTRUE
  Bool_t bSetSameChargePositiv;   	// used if bDifferentCharge: if kTRUE use positiv, if kFALSE use negative (default kTRUE)
  Int_t fMixedHarmonic;			// Harmonic of special mixed particle analysis
  TH1F *fCounterHistogram;       	//! for some checks
  TProfile *fProfileEventCuts;  	//! Profile to save the cut values for event selection
  TProfile *fProfileTrackCuts;  	//! Profile to save the cut values for track selection
  TList *fFinalResultsList[16];      	//! List to hold all histograms with final results for a specific centrality bin. Up to 16 centraliy bins possible

  

  ClassDef(AliAnalysisTaskStudentsML,39); 

};

//================================================================================================================

#endif
