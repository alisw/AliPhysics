/*
  Example of AliEffFDContainer usage.
  The PCC framework (AliMCSpectraWeights) used in AliEffFDContainer is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/
#ifndef ALIANALYSISTASKEFFANDFDTEST__H
#define ALIANALYSISTASKEFFANDFDTEST__H
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "TString.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"


class TList;
class AliVTrack;
class AliVVertex;
class AliInputEventHandler;
class AliAnalysisUtils;
class AliVParticle;
class AliEffFDContainer;

class AliAnalysisTaskEffFDExample : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEffFDExample();
  AliAnalysisTaskEffFDExample(const char *name, Bool_t IsMC=kTRUE, TString pf="");
  virtual ~AliAnalysisTaskEffFDExample();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void SetEta(Double_t etaMin, Double_t etaMax) { fEtaMin=etaMin; fEtaMax = etaMax; };
  void SetContPF(TString newval) { fContPF = newval; };
  void SetPtBins(Int_t nBins, Double_t *ptbins) { if(fPtAxis) delete fPtAxis; fPtAxis = new TAxis(nBins,ptbins); };
  void SetMultiBins(Int_t nBins, Double_t *multibins) { if(fMultiAxis) delete fMultiAxis; fMultiAxis = new TAxis(nBins,multibins); };
  void SetCentEstimator(TString newval) {fCentEst = newval; };
  void ClearTCList() { l_ClearTCList(); fFBtoAdd=0; };
  void AddTrackCut(AliESDtrackCuts* incut) { l_CreateTCList(); fTCtoAdd->Add(incut); };
  void AddTrackCut(Int_t fb) { fFBtoAdd+=fb; };
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskEffFDExample(const AliAnalysisTaskEffFDExample&);
  AliAnalysisTaskEffFDExample& operator=(const AliAnalysisTaskEffFDExample&);
  Bool_t fIsMC;
  AliMCEvent *fMCEvent; //! MC event
  UInt_t fTriggerType;
  AliEffFDContainer *fEfFd;
  Double_t fEtaMin;
  Double_t fEtaMax;
  TString fContPF;
  TAxis *fPtAxis; //for storing pT axis
  TAxis *fMultiAxis; //for storing cent/multi axis
  Double_t *fPtBins; //!
  Int_t fNPtBins; //!
  Double_t *fMultiBins; //!
  Int_t fNMultiBins; //!
  TString fCentEst;
  TList *fTCtoAdd;
  Int_t fFBtoAdd;
  Double_t *GetBinsFromAxis(TAxis *inax);
  void l_CreateTCList() { if(fTCtoAdd) return; fTCtoAdd = new TList(); fTCtoAdd->SetOwner(kTRUE); };
  void l_ClearTCList() { if(fTCtoAdd) fTCtoAdd->Clear(); };
  ClassDef(AliAnalysisTaskEffFDExample,1);
};

#endif
