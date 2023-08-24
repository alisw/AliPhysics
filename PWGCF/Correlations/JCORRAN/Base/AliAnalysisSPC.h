/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/********************************************** 
* template class for student projects         *
* author: Marcel Lesch (marcel.lesch@cern.ch) *
**********************************************/ 

#ifndef ALIANALYSISSPC_H
#define ALIANALYSISSPC_H

#include "TProfile.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TComplex.h"
#include <TArrayD.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include "AliHeader.h"
#include "TExMap.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TFile.h"
#include "AliJEfficiency.h"
#include "TSystem.h"

class TClonesArray;

//================================================================================================================

class AliAnalysisSPC{
public:
  AliAnalysisSPC();
  AliAnalysisSPC(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisSPC(); 
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  // 0.) Methods called in the constructor:
  virtual void InitializeArrays();
 
  // 1.) Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();

  // 2.) Methods called in UserExec(Option_t *) or subsquent Methods:
  virtual void PhysicsAnalysis();

  //a) Methods used to assure global quality and track selection
  Int_t SelectCentrality(Double_t CentralityValue);
  virtual void FisherYatesRandomizing(Int_t Mult, Int_t *RandomIndex); 

  //b) Methods used for Multi-Particle Correlation Techniques
  virtual void CalculateQvectors(Int_t CalculateQvectors_nParticles, Double_t* CalculateQvectors_angles, Double_t* CalculateQvectors_weights);
  TComplex CalculateMixedQVectors(Double_t Harm, Int_t M_A, Int_t M_B, Double_t* Ang_A, Double_t* Ang_B);
  TComplex Q(Int_t n, Int_t p);
  TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip);
  virtual void Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8, Int_t h9, Int_t h10, Int_t h11, Int_t h12, Int_t h13, Int_t h14, Double_t *Correlation_Data);

  virtual void MainTask(Int_t MainTask_CentBin, Int_t MainTask_Mult, Double_t* MainTask_Angle_Array, Double_t* MainTask_Weight_Array);
  virtual void MixedParticle(Int_t MP_CentBin, Int_t Harmonicus, Int_t Mixed_Mult_A, Double_t* Mixed_Angle_A, Int_t Mixed_Mult_B, Double_t* Mixed_Angle_B);
 
  virtual void ComputeTPCWithEtaGaps(Int_t CentralityBin, Int_t numberOfParticles, Double_t* angles, Double_t* pWeights, Double_t* pseudorapidity);

  // 3.) Methods called in Terminate():
  // ...

  // 4.) Setters and getters: 
  void SetInputList(TClonesArray *inputarray){fInputList = inputarray;}
  TClonesArray *GetInputList() const{return fInputList;}

  void SetEventCentrality( float cent ){fCentrality = cent;}

  TList* GetMainList() const{return fHistList;} // get the list for external task

  void SetDebugLevel(int debuglevel){fDebugLevel = debuglevel; cout <<"setting Debug Level = " << fDebugLevel << endl;}

  void SetSaveAllQA(Bool_t SaveQA){this->bSaveAllQA=SaveQA;}

  void SetUseWeights(Bool_t WeightsNUE, Bool_t WeightsNUA){this->bUseWeightsNUE = WeightsNUE; this->bUseWeightsNUA = WeightsNUA;}

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
	this->bSetSameChargePositive = PositiveCharge;	}

  void SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4, Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10, Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16 )
 {this->fcent_0 = cen0; this->fcent_1 = cen1; this->fcent_2 = cen2; this->fcent_3 = cen3; this->fcent_4 = cen4; this->fcent_5 = cen5; this->fcent_6 = cen6; this->fcent_7 = cen7; this->fcent_8 = cen8; this->fcent_9 = cen9; this->fcent_10 = cen10; this->fcent_11 = cen11; this->fcent_12 = cen12; this->fcent_13 = cen13; this->fcent_14 = cen14; this->fcent_15 = cen15; this->fcent_16 = cen16;} 

  void SetInitializeCentralityArray(); //Set Centrality array 

  void SetEtaGaps(Bool_t ComputeEtaGap, Float_t EtaGap)
  {this->bComputeEtaGap = ComputeEtaGap; this->fEtaGap = EtaGap; } 

private:
  AliAnalysisSPC(const AliAnalysisSPC& aatmpf);
  AliAnalysisSPC& operator=(const AliAnalysisSPC& aatmpf);
  
  // 0.) Base lists:
  TList *fHistList; // base list to hold all output object (a.k.a. grandmother of all lists)

  TClonesArray *fInputList;
  Double_t fCentrality;
  Int_t fDebugLevel;

  // 1.) Control histograms: 
  TList *fCentralityList[16];		//! Will be one list per certain centrality bin. Up to 16 centraliy bins possible
  TList *fControlHistogramsList[16]; 	//! List to hold all control histograms for a specific centrality bin. Up to 16 centraliy bins possible
 
  TH1F *fPTHistogram[16][2]; 		//! 0: P_t After Track Selection, 1: P_t After Track Selection (second)
  TH1F *fPhiHistogram[16][2]; 		//! 0: Phi After Track Selection, 1: Phi After Track Selection (second)
  TH1F *fEtaHistogram[16][2]; 		//! 0: Eta After Track Selection, 1: Eta After Track Selection (second)
  TH1F *fMultHistogram[16][2]; 		//! 0: Mult. After Track Selection 1: Mult. After Track Selection (second)
  TH1I *fChargeHistogram[16]; 	//! Charge After Track Selection 
  TH1F *fCentralityHistogram[16]; 	//! Centrality After Corresponding Cut
  TProfile *fPhiWeightProfile[16];	//! Profile to check phi weights

  //2.) SelectionCuts
  Bool_t bSaveAllQA;			// if kTRUE: All Standard QA Histograms are saved (default kTRUE)

  //Centrality
  Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
  					//fcent_i holds the edge of a centrality bin
  Int_t fCentralityBins;		//will be set to 16, for at maximum 16 bins in case of centrality 0 to 80 in 5% steps. Less bins and different steps may be used
  Float_t fcentralityArray[17];		//will hold our edges of the centrality bins

  //Fisher-Yales
  Bool_t bDoFisherYates;		//if kTRUE: Do Fisher Yates Mixing of phi, pt and eta arrays after track selection (default: kFALSE)
  Float_t fFisherYatesCutOff;		//How much percentage of the orginal particles are kept, e.g. if 0.7 only 70% of the current particles are kept for analysis

  //Weights
  Bool_t bUseWeightsNUE; 
  Bool_t bUseWeightsNUA; 

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
  TProfile *fResultsAlternativeError[16]; //! final centrality result
  TProfile *fCovResults[16];         	//! TProfile to store terms needed for Covariance 
  TProfile *fJoinedCovResults[16];      //! TProfile to store joined Covariance term calculated as one correlator <z> instead of product of two correlators <x*y> 
  TProfile *fMixedParticleHarmonics[16];//! Stores output for special mixed particle analysis
  Bool_t bDoMixed;		 	// if kTRUE: Do special mixed particle analysis, default kFALSE (MainTask)
  Bool_t bDifferentCharge; 	 	// used in DoMixed: if kTRUE mixed particle analysis between positiv and negativ
				 	//		    if kFALSE mixed particle analysis between same charge 
				 	//		    (only positiv or only negativ particles)
				 	// Default kTRUE
  Bool_t bComputeEtaGap;		// Do eta gap computation if kTRUE. Default kFALSE
  Float_t fEtaGap;			// Value of eta gap
  TProfile *fProfileTPCEta[16];		//! Profile for 2-particle eta gap computation
  

  Bool_t bSetSameChargePositive;   	// used if bDifferentCharge: if kTRUE use positiv, if kFALSE use negative (default kTRUE)
  Int_t fMixedHarmonic;			// Harmonic of special mixed particle analysis
  TH1F *fCounterHistogram;       	//! for some checks
  TProfile *fProfileTrackCuts;  	//! Profile to save the cut values for track selection
  TList *fFinalResultsList[16];      	//! List to hold all histograms with final results for a specific centrality bin. Up to 16 centraliy bins possible

  ClassDef(AliAnalysisSPC,5); 
};

//================================================================================================================

#endif
