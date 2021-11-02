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
#include "AliAnalysisTaskHOCFA.h"
#include "TDirectory.h"
#include "TComplex.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliJCatalystTask.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

class TClonesArray;
class AliJFlowHistos;
using namespace std;

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

// Methods specific to this class.
  void BookHistos(TClonesArray *inList);  // TBI: There is no body for this method?
  bool IsMC() const {return fIsMC;}
  void SetIsMC(bool b) {fIsMC = b;}
  AliJCatalystTask *GetJCatalystTask() {return fJCatalystTask;}
  void SetJCatalystTaskName(TString name) {fJCatalystTaskName = name;}
  TString GetJCatalystTaskName() {return fJCatalystTaskName;}

// Methods from the analysis task.
  void HOCFASetDebugLevel(int debug) {fHOCFADebugLevel = debug;}
  void HOCFASetCentralityBinning(int nBins) {fHOCFAnCentralityBins = nBins;}
  void HOCFASetCentralityArray(TString values) {fHOCFAvalues = values;}
  void HOCFASetMinMultiplicity(int minMult) {fHOCFAMultiplicityMin = minMult;}
  void HOCFASetPtRange(double minPt, double maxPt) {
    fHOCFAPtMin = minPt; fHOCFAPtMax = maxPt;
  }
  void HOCFASetParticleWeights(bool weightsNUE, bool weightsNUA) {
    fHOCFAUseWeightsNUE = weightsNUE; fHOCFAUseWeightsNUA = weightsNUA;
  }
  void HOCFASetObservable(bool observ) {fHOCFAGetSC3h = observ;}
  void HOCFASetNumberCombi(int combi) {fHOCFANCombi = combi;}
  void HOCFASetHarmoArray(TString combiString) {fHOCFAcombi = combiString;}

private:
  AliJCatalystTask *fJCatalystTask;   // Pointer to the catalyst task.
  TString fJCatalystTaskName;         // Name of the catalyst task.
  bool fIsMC;                         // MC or real data.
  AliAnalysisTaskHOCFA *fHOCFATask;   // Pointer to the analysis task.

  int fHOCFADebugLevel;               // Select how much is printed in the terminal.
  int fHOCFAnCentralityBins;          //! Number of centrality bins (Size(array)-1).
  int fHOCFAMultiplicityMin;          // Minimum multiplicity to have valid events.
  double fHOCFAPtMin;                 // Minimum transverse momentum.
  double fHOCFAPtMax;                 // Maximum transverse momentum.
  bool fHOCFAUseWeightsNUE;           // kTRUE: Enable the non-unit NUE corrections.
  bool fHOCFAUseWeightsNUA;           // kTRUE: Enable the non-unit NUA corrections.
  bool fHOCFAGetSC3h;                 // kTRUE: Calculate SC(k,l,m), else AC(m,n).
  int fHOCFANCombi;                   // Number of combinations of harmonics (max 6).
  TString fHOCFAvalues;               // Values for the centrality edges (max 17).
  TString fHOCFAcombi;                // Values for the harmonics combinations.

  ClassDef(AliJHOCFATask, 3);
};

#endif