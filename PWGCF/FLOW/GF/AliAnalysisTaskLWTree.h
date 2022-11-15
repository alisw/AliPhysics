/*
  Example of AliEffFDContainer usage.
  The PCC framework (AliMCSpectraWeights) used in AliEffFDContainer is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/
#ifndef ALIANALYSISTASKLWTREE__H
#define ALIANALYSISTASKLWTREE__H
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliGFWFilter.h"
#include "AliLWUtils.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "AliAODForwardMult.h"

class TList;
class AliInputEventHandler;
class AliAnalysisUtils;
class AliLWTrack;
class AliLWEvent;

class AliAnalysisTaskLWTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLWTree();
  AliAnalysisTaskLWTree(const char *name);
  virtual ~AliAnalysisTaskLWTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void SetEta(Double_t etaMin, Double_t etaMax) { fEtaMin=etaMin; fEtaMax = etaMax; };
  void SetCentEstimator(TString newval) {fCentEst = newval; };
  void ClearTCList() { l_ClearTCList(); fFBtoAdd=0; };
  void AddTrackCut(AliESDtrackCuts* incut) { l_CreateTCList(); fTCtoAdd->Add(incut); };
  void AddTrackCut(Int_t fb) { fFBtoAdd+=fb; };
  void SetVertexCut(Double_t newval) {fVtxZCut = newval; };
  void SetupFlags();
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskLWTree(const AliAnalysisTaskLWTree&);
  AliAnalysisTaskLWTree& operator=(const AliAnalysisTaskLWTree&);
  TTree *fOutTree;
  AliLWEvent *fLWEvent;
  TClonesArray *fTPCTracks;
  TClonesArray *fFMDTracks;
  UInt_t fTriggerType;
  Double_t fVtxZCut;
  Double_t fEtaMin;
  Double_t fEtaMax;
  TString fCentEst;
  TList *fTCtoAdd;
  Int_t fFBtoAdd;
  UInt_t fEvNomFlag; //Relevant for AODs
  UInt_t fTrNomFlag; //Relevant for AODs
  void l_CreateTCList() { if(fTCtoAdd) return; fTCtoAdd = new TList(); fTCtoAdd->SetOwner(kTRUE); };
  void l_ClearTCList() { if(fTCtoAdd) fTCtoAdd->Clear(); };
  TRandom gRndm;
  ClassDef(AliAnalysisTaskLWTree,1);
};

#endif
