#ifndef ALIJSPCTASK_H
#define ALIJSPCTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations 
// author: D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla 
// Finland 
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <TDirectory.h>
#include <TComplex.h>
#include <AliLog.h>
#include <AliAnalysisTaskSE.h>
#include <AliJCatalystTask.h>
#include "AliJEfficiency.h"
#include "AliAnalysisSPC.h"


using namespace std;
//==============================================================
class TClonesArray;
class AliJFlowHistos;

class AliJSPCTask : public AliAnalysisTaskSE {

 public:
  AliJSPCTask();
  AliJSPCTask(const char *name);
  AliJSPCTask(const AliJSPCTask& ap);   
  AliJSPCTask& operator = (const AliJSPCTask& ap);
  virtual ~AliJSPCTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  void BookHistos(TClonesArray *inList);
  Bool_t IsMC()const{ return fIsMC; }
  void SetIsMC(Bool_t b){ fIsMC=b; }

  //Task Specific Setters
  void JSPCSetSaveAllQA(Bool_t SaveQA){bJSPCSaveAllQA=SaveQA;}
  void JSPCSetUseWeights(Bool_t WeightsNUE, Bool_t WeightsNUA){bJSPCUseWeightsNUE = WeightsNUE; bJSPCUseWeightsNUA = WeightsNUA; }
  void JSPCSetFisherYates(Bool_t DoFY, Float_t CutOff)
  { bJSPCDoFisherYates=DoFY; fJSPCFisherYatesCutOff=CutOff; } 
  void JSPCSetMinNuPar(Int_t top){fJSPCMinNumberPart = top;} 
  Int_t JSPCGetMinNuPar() const {return fJSPCMinNumberPart;}

   void JSPCSetCorrSet1(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {fJSPCNumber=Number; fJSPCa1=a; fJSPCa2=b; fJSPCa3=c; fJSPCa4=d; fJSPCa5=e; fJSPCa6=f; fJSPCa7=g;}
  void JSPCSetCorrSet2( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { fJSPCNumberSecond=Number; fJSPCb1=a; fJSPCb2=b; fJSPCb3=c; fJSPCb4=d; fJSPCb5=e; fJSPCb6=f; fJSPCb7=g;}
  void JSPCSetCorrSet3( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { fJSPCNumberThird=Number; fJSPCd1=a; fJSPCd2=b; fJSPCd3=c; fJSPCd4=d; fJSPCd5=e; fJSPCd6=f; fJSPCd7=g;}
  void JSPCSetCorrSet4( Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  { fJSPCNumberFourth=Number; fJSPCe1=a; fJSPCe2=b; fJSPCe3=c; fJSPCe4=d; fJSPCe5=e; fJSPCe6=f; fJSPCe7=g;}
  void JSPCSetCorrSet5(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {fJSPCNumberFifth=Number; fJSPCf1=a; fJSPCf2=b; fJSPCf3=c; fJSPCf4=d; fJSPCf5=e; fJSPCf6=f; fJSPCf7=g;}
  void JSPCSetCorrSet6(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {fJSPCNumberSixth=Number; fJSPCg1=a; fJSPCg2=b; fJSPCg3=c; fJSPCg4=d; fJSPCg5=e; fJSPCg6=f; fJSPCg7=g;}
  void JSPCSetCorrSet7(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {fJSPCNumberSeventh=Number; fJSPCh1=a; fJSPCh2=b; fJSPCh3=c; fJSPCh4=d; fJSPCh5=e; fJSPCh6=f; fJSPCh7=g;}
  void JSPCSetCorrSet8(Int_t Number, Int_t a, Int_t b, Int_t c, Int_t d, Int_t e, Int_t f, Int_t g)
  {fJSPCNumberEighth=Number; fJSPCi1=a; fJSPCi2=b; fJSPCi3=c; fJSPCi4=d; fJSPCi5=e; fJSPCi6=f; fJSPCi7=g;}
  void JSPCSetMixed(Bool_t top, Int_t nop, Bool_t DifferentCharge, Bool_t PositiveCharge)
  {bJSPCDoMixed = top; fJSPCMixedHarmonic = nop; bJSPCDifferentCharge = DifferentCharge; bJSPCSetSameChargePositive = PositiveCharge;}

  void JSPCSetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4, Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10, Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16 )
 {fJSPCcent_0 = cen0; fJSPCcent_1 = cen1; fJSPCcent_2 = cen2; fJSPCcent_3 = cen3; fJSPCcent_4 = cen4; fJSPCcent_5 = cen5; fJSPCcent_6 = cen6; fJSPCcent_7 = cen7; fJSPCcent_8 = cen8; fJSPCcent_9 = cen9; fJSPCcent_10 = cen10; fJSPCcent_11 = cen11; fJSPCcent_12 = cen12; fJSPCcent_13 = cen13; fJSPCcent_14 = cen14; fJSPCcent_15 = cen15; fJSPCcent_16 = cen16;} 
  void JSPCSetInitializeCentralityArray(); //Set Centrality array 

  void JSPCSetEtaGaps(Bool_t ComputeEtaGap, Float_t EtaGap)
  {this->bJSPCComputeEtaGap = ComputeEtaGap; this->fJSPCEtaGap = EtaGap; } 

  // Methods specific for this class
  void SetJCatalystTaskName(TString name){ fJCatalystTaskName=name; } // Setter for filter task name
  TString GetJCatalystTaskName(){ return fJCatalystTaskName; } // Setter for filter task name

 private:

  AliJCatalystTask *fJCatalystTask;  //
  TString           fJCatalystTaskName; // Name for JCatalyst task
  Bool_t      fIsMC;       // MC data or real data
  AliAnalysisSPC *fSPC;

  //Saving Minor QA
  Bool_t bJSPCSaveAllQA;      // if kTRUE: All Standard QA Histograms are saved (default kTRUE)

  //Centrality
  Float_t fJSPCcent_0, fJSPCcent_1, fJSPCcent_2, fJSPCcent_3, fJSPCcent_4, fJSPCcent_5, fJSPCcent_6, fJSPCcent_7, fJSPCcent_8, fJSPCcent_9, fJSPCcent_10, fJSPCcent_11, fJSPCcent_12, fJSPCcent_13, fJSPCcent_14, fJSPCcent_15, fJSPCcent_16; //fcent_i holds the edge of a centrality bin
 
  //Fisher-Yates
  Bool_t bJSPCDoFisherYates;    //if kTRUE: Do Fisher Yates Mixing of phi, pt and eta arrays after track selection (default: kFALSE)
  Float_t fJSPCFisherYatesCutOff;   //How much percentage of the orginal particles are kept, e.g. if 0.7 only 70% of the current particles are kept for analysis

  //Minimum Number of Particles for Correlation
  Int_t fJSPCMinNumberPart;

  //Weights
  Bool_t bJSPCUseWeightsNUE; 
  Bool_t bJSPCUseWeightsNUA;
 
  //Correlators
  Int_t fJSPCNumber;              // Number of correlation first correlator
  Int_t fJSPCNumberSecond;            // Number of correlation second correlator
  Int_t fJSPCNumberThird;             // Number of correlation third correlator
  Int_t fJSPCNumberFourth;      // Number of correlation fourth correlator
  Int_t fJSPCNumberFifth;     // Number of correlation fifth correlator
  Int_t fJSPCNumberSixth;     // Number of correlation sixth correlator
  Int_t fJSPCNumberSeventh;     // Number of correlation seventh correlator
  Int_t fJSPCNumberEighth;      // Number of correlation eigth correlator

  Int_t fJSPCa1, fJSPCa2, fJSPCa3, fJSPCa4, fJSPCa5, fJSPCa6, fJSPCa7; //first set of harmonics
  Int_t fJSPCb1, fJSPCb2, fJSPCb3, fJSPCb4, fJSPCb5, fJSPCb6, fJSPCb7; //second set of harmonics
  Int_t fJSPCd1, fJSPCd2, fJSPCd3, fJSPCd4, fJSPCd5, fJSPCd6, fJSPCd7; //third set of harmonics
  Int_t fJSPCe1, fJSPCe2, fJSPCe3, fJSPCe4, fJSPCe5, fJSPCe6, fJSPCe7; //fourth set of harmonics
  Int_t fJSPCf1, fJSPCf2, fJSPCf3, fJSPCf4, fJSPCf5, fJSPCf6, fJSPCf7;  //sixth set of harmonics 
  Int_t fJSPCg1, fJSPCg2, fJSPCg3, fJSPCg4, fJSPCg5, fJSPCg6, fJSPCg7;  //seventh set of harmonics 
  Int_t fJSPCh1, fJSPCh2, fJSPCh3, fJSPCh4, fJSPCh5, fJSPCh6, fJSPCh7;  //harmonics
  Int_t fJSPCi1, fJSPCi2, fJSPCi3, fJSPCi4, fJSPCi5, fJSPCi6, fJSPCi7;  //eigth set of harmonics 

  //Mixed Charge
  Bool_t bJSPCDoMixed;      // if kTRUE: Do special mixed particle analysis, default kFALSE (MainTask)
  Bool_t bJSPCDifferentCharge;    // used in DoMixed: if kTRUE mixed particle analysis between positiv and negativ
          //        if kFALSE mixed particle analysis between same charge 
          //        (only positiv or only negativ particles)
          // Default kTRUE
  Bool_t bJSPCSetSameChargePositive;     // used if bDifferentCharge: if kTRUE use positiv, if kFALSE use negative (default kTRUE)
  Int_t fJSPCMixedHarmonic;     // Harmonic of special mixed particle analysis

  Bool_t bJSPCComputeEtaGap;		// Do eta gap computation if kTRUE. Default kFALSE
  Float_t fJSPCEtaGap;			// Value of eta gap

  ClassDef(AliJSPCTask, 3); 
};
#endif // AliJSPCTask_H
