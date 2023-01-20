/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Interface task to get the tracks from the catalyst to the analysis task.    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ -------------------------------------------------------------------------- */
#ifndef AliJHOCFATask_H
#define AliJHOCFATask_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "AliAnalysisTaskHOCFA.h"
#include "AliAnalysisTaskSE.h"
#include "AliJCatalystTask.h"
#include "AliLog.h"
#include "TComplex.h"
#include "TDirectory.h"

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
  void HOCFASetCentralityArray(TString values) {fHOCFAvalues = values;}
  void HOCFASetCentralityBinning(int nBins) {fHOCFAnCentralityBins = nBins;}
  void HOCFASetMinMultiplicity(int minMult) {fHOCFAMultiplicityMin = minMult;}

  void HOCFASetPtRange(double minPt, double maxPt) {fHOCFAPtMin = minPt; fHOCFAPtMax = maxPt;}  
  void HOCFASetEtaGap(bool etaGap, float myGap) {fHOCFAApplyEtaGap = etaGap; fHOCFAEtaGap = myGap;}
  void HOCFASetParticleWeights(bool weightsNUE, bool weightsNUA) {
    fHOCFAUseWeightsNUE = weightsNUE; fHOCFAUseWeightsNUA = weightsNUA;
  }
  void HOCFASetCentralityWeights(bool weightsCent) {fHOCFAUseWeightsCent = weightsCent;}
  void HOCFASetObservable(bool myObs, bool myOrder) {fHOCFAGetSC = myObs; fHOCFAGetLower = myOrder;}

 private:
  AliJCatalystTask *fJCatalystTask;   // Pointer to the catalyst task.
  TString fJCatalystTaskName;         // Name of the catalyst task.
  AliAnalysisTaskHOCFA *fHOCFATask;   // Pointer to the analysis task.
  bool fIsMC;                         // MC or real data.

  int fHOCFADebugLevel;               // Select how much is printed in the terminal.

  TString fHOCFAvalues;               // String gathering all the values for the centrality binning.
  int fHOCFAnCentralityBins;          //! Number of centrality bins (Size(array)-1).
  int fHOCFAMultiplicityMin;          // Minimum multiplicity to have valid events.

  double fHOCFAPtMin;                 // Minimum transverse momentum.
  double fHOCFAPtMax;                 // Maximum transverse momentum.
  float fHOCFAEtaGap;                 // Value of the gap (default: 0.).
  bool fHOCFAApplyEtaGap;             // kTRUE: Get the 2p correlators with an eta gap.
  bool fHOCFAUseWeightsNUE;           // kTRUE: Enable the non-unit NUE corrections.
  bool fHOCFAUseWeightsNUA;           // kTRUE: Enable the non-unit NUA corrections.
  bool fHOCFAUseWeightsCent;          // kTRUE: Enable the non-unit centrality corrections for LHC15o.

  bool fHOCFAGetSC;                   // kTRUE: Measure 2-h and 3-h SC, else 2-h AC.
  bool fHOCFAGetLower;                // kTRUE: Measure the terms for the lower harmonics.

  ClassDef(AliJHOCFATask, 6);
};

#endif  // AliJHOCFATask_H