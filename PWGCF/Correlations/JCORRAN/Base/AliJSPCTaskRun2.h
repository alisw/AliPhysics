/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the numerator and denominators for the /
/ various SPC measured for Run2 Pb-Pb data                                    /
/                                                                             /
/ Authors: Cindy Mordasini (cindy.mordasini@cern.ch)                          /
/          Maxim Virta                                                        /
/ -------------------------------------------------------------------------- */
#ifndef ALIJSPCTASKRUN2_H
#define ALIJSPCTASKRUN2_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TComplex.h>
#include <TDirectory.h>

#include "AliAnalysisSPCRun2.h"
#include <AliAnalysisTaskSE.h>
#include <AliJCatalystTask.h>
#include "AliJEfficiency.h"
#include <AliLog.h>

using namespace std;

class TClonesArray;
class AliJFlowHistos;

class AliJSPCTaskRun2 : public AliAnalysisTaskSE {
 public:
  AliJSPCTaskRun2();
  AliJSPCTaskRun2(const char *name);
  AliJSPCTaskRun2(const AliJSPCTaskRun2& ap);   
  AliJSPCTaskRun2& operator = (const AliJSPCTaskRun2& ap);
  virtual ~AliJSPCTaskRun2();

  // Methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* option);

  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  void BookHistos(TClonesArray *inList);
  Bool_t IsMC() const {return fIsMC;}
  void SetIsMC(Bool_t b) {fIsMC=b;}

  // Task specific setters
  void AliSPCRun2SetSaveAllQA(Bool_t SaveQA) {bJSPCRun2SaveAllQA=SaveQA;}
  void AliSPCRun2SetUseWeights(Bool_t WeightsNUE, Bool_t WeightsNUA) {
    bAliSPCRun2UseWeightsNUE = WeightsNUE;
    bAliSPCRun2UseWeightsNUA = WeightsNUA;
  }
  void AliSPCRun2SetMinNuPar(Int_t top) {fAliSPCRun2MinNumberPart = top;} 
  Int_t AliSPCRun2GetMinNuPar() const {return fAliSPCRun2MinNumberPart;}

  void AliSPCRun2SetCorrSet(Int_t obsInd, Int_t harmo[8]) {
    for (int i = 0; i < 8; i++) {fAliSPCRun2HarmosArray[obsInd][i] = harmo[i];}
  }

  void AliSPCRun2SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3,
    Float_t cen4, Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9,
    Float_t cen10, Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15,
    Float_t cen16) {
    fAliSPCRun2cent_0 = cen0; fAliSPCRun2cent_1 = cen1; fAliSPCRun2cent_2 = cen2;
    fAliSPCRun2cent_3 = cen3; fAliSPCRun2cent_4 = cen4; fAliSPCRun2cent_5 = cen5;
    fAliSPCRun2cent_6 = cen6; fAliSPCRun2cent_7 = cen7; fAliSPCRun2cent_8 = cen8;
    fAliSPCRun2cent_9 = cen9; fAliSPCRun2cent_10 = cen10; fAliSPCRun2cent_11 = cen11;
    fAliSPCRun2cent_12 = cen12; fAliSPCRun2cent_13 = cen13; fAliSPCRun2cent_14 = cen14;
    fAliSPCRun2cent_15 = cen15; fAliSPCRun2cent_16 = cen16;} 
  void AliSPCRun2SetInitializeCentralityArray();

  void AliSPCRun2SetEtaGaps(Bool_t ComputeEtaGap, Float_t EtaGap) {
    bAliSPCRun2ComputeEtaGap = ComputeEtaGap; fAliSPCRun2EtaGap = EtaGap;} 

  void SetJCatalystTaskName(TString name) {fJCatalystTaskName = name;}
  TString GetJCatalystTaskName() {return fJCatalystTaskName;}

 private:
  AliJCatalystTask *fJCatalystTask;   // Instance of the catalyst task.
  TString fJCatalystTaskName;         // Name for JCatalyst task
  Bool_t fIsMC;                       // MC data or real data
  AliAnalysisSPCRun2 *fSPC;           // Instance of the analysis task.

  // Saving Minor QA.
  Bool_t bAliSPCRun2SaveAllQA;        // kTRUE: Save the standard QA histograms (default: kTRUE).

  // Centrality.
  Float_t fAliSPCRun2cent_0, fAliSPCRun2cent_1, fAliSPCRun2cent_2, fAliSPCRun2cent_3,
    fAliSPCRun2cent_4, fAliSPCRun2cent_5, fAliSPCRun2cent_6, fAliSPCRun2cent_7, fAliSPCRun2cent_8,
    fAliSPCRun2cent_9, fAliSPCRun2cent_10, fAliSPCRun2cent_11, fAliSPCRun2cent_12,
    fAliSPCRun2cent_13, fAliSPCRun2cent_14, fAliSPCRun2cent_15, fAliSPCRun2cent_16;

  Int_t fAliSPCRun2MinNumberPart;     // Minimum Number of Particles for Correlation.

  Bool_t bAliSPCRun2UseWeightsNUE;    // kTRUE: Use non-unit NUE weights.
  Bool_t bAliSPCRun2UseWeightsNUA;    // kTRUE: Use non-unit NUA weights.
 
  Bool_t bJSPCRun2SaveAllQA;
  Int_t fAliSPCRun2HarmosArray[12][8];  // Array of combinations of harmonics for the SPC.
    // Can deal with maximum 10 different SPC.

  Bool_t bAliSPCRun2ComputeEtaGap;    // Do eta gap computation if kTRUE. Default kFALSE
  Float_t fAliSPCRun2EtaGap;          // Value of eta gap

  ClassDef(AliJSPCTaskRun2, 1); 
};

#endif  // ALIJSPCTASKRUN2_H
