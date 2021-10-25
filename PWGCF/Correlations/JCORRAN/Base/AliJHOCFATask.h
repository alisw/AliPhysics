#ifndef AliJHOCFATask_H
#define AliJHOCFATask_H

/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Interface task to get the tracks from the catalyst to the analysis task.    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ -------------------------------------------------------------------------- */
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
#include "AliAnalysisTaskHOCFA.h"

using namespace std;
class TClonesArray;
class AliJFlowHistos;

class AliJHOCFATask : public AliAnalysisTaskSE {
public:
// Methods inherited from AliAnalysisTaskSE.
  AliJHOCFATask();
  AliJHOCFATask(const char *name);
  AliJHOCFATask(const AliJHOCFATask& ap);   
  AliJHOCFATask& operator = (const AliJHOCFATask& ap);
  virtual ~AliJHOCFATask();

  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* );

// Methods specific for this class.
  void BookHistos(TClonesArray *inList);  // TBI: There is no body for this method?
  Bool_t IsMC() const {return fIsMC;}
  void SetIsMC(Bool_t b) {fIsMC = b;}
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  void SetJCatalystTaskName(TString name) {fJCatalystTaskName = name;}
  TString GetJCatalystTaskName() {return fJCatalystTaskName;}

// Methods from the analysis task.
  void HOCFASetCentralityBinning(Int_t nBins) {fHOCFANCentralityBins = nBins;}
  void HOCFASetCentralityArray(TString values) {fHOCFAvalues = values;}
  void HOCFASetMinMultiplicity(Int_t minMult) {fHOCFAMultiplicityMin = minMult;}
  void HOCFASetParticleWeights(Bool_t weightsNUE, Bool_t weightsNUA) {fHOCFAUseWeightsNUE = weightsNUE; fHOCFAUseWeightsNUA = weightsNUA;}
  void HOCFASetNumberCombi(Int_t combi) {fHOCFANCombi = combi;}
  void HOCFASetHarmoArray(TString combiString){fHOCFAcombi = combiString;}

private:
  AliJCatalystTask *fJCatalystTask; // Pointer to the catalyst task.
  TString fJCatalystTaskName; // Name of the catalyst task.
  Bool_t fIsMC; // MC or real data.
  AliAnalysisTaskHOCFA *fHOCFATask; // Pointer to the analysis task.

  Int_t fHOCFANCentralityBins; //! Number of centrality bins in the division (Size(array)-1).
  Int_t fHOCFAMultiplicityMin; // Minimum multiplicity to calculate the correlators.
  Bool_t fHOCFAUseWeightsNUE; // kTrue: Enable the non-unit NUE corrections.
  Bool_t fHOCFAUseWeightsNUA; // kTrue: Enable the non-unit NUA corrections.
  Int_t fHOCFANCombi;  // Number of combinations of harmonics (max 6).
  TString fHOCFAvalues; // Values for the centrality edges (max 17).
  TString fHOCFAcombi;  // Values for the harmonics combinations.

  ClassDef(AliJHOCFATask, 2);
};

#endif
