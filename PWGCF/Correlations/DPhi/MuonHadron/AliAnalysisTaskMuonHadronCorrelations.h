#ifndef AliAnalysisTaskMuonHadronCorrelations_H
#define AliAnalysisTaskMuonHadronCorrelations_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TH2D.h"

//====================================================================================================================================================

class  AliAnalysisTaskMuonHadronCorrelations : public AliAnalysisTaskSE {

 public:
 
  enum {kSingleEvent, kMixedEvent};

  AliAnalysisTaskMuonHadronCorrelations();
  AliAnalysisTaskMuonHadronCorrelations(const Char_t *name);
  virtual ~AliAnalysisTaskMuonHadronCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // ------------- Cuts -----------------

  void SetFilterBitCentralBarrel(Int_t filter) { fFilterBitCentralBarrel = filter; }
  void SetMaxEtaCentralBarrel(Double_t eta) { fMaxEtaCentralBarrel = eta; }
  void SetTriggerMatchLevelMuon(Short_t level) { fTriggerMatchLevelMuon = level; }
  //  void SetMaxChi2Muon(Double_t chi2Max) { fMaxChi2Muon = chi2Max; }
  void SetRAbsRangeMuon (Double_t rAbsMin,Double_t rAbsMax) { fMinRAbsMuon = rAbsMin; fMaxRAbsMuon = rAbsMax; }

  void SetTriggerWord(TString triggerWord) { fTriggerWord = triggerWord; fIsTriggerSet = kTRUE; }

  // ------------- Analysis -------------

  Float_t GetV0Multiplicity();
  Double_t GetITSMultiplicity();
  Bool_t IsTriggerFired();
  TClonesArray* GetAcceptedTracksCentralBarrel(AliAODEvent *aodEvent);
  TClonesArray* GetAcceptedTracksMuonArm(AliAODEvent *aodEvent);
  void SetPtBinning(Int_t nBins, Double_t *limits);
  void SetCentBinning(Int_t nBins, Double_t *limits);
  void SetMultBinning(Int_t nBins, Double_t *limits);
  void SetCentMethod(const Char_t *method) { fCentMethod = method; }
  void FillHistograms(Int_t centrality, Int_t option);
  Int_t GetCentBin();
  Int_t GetMultBin();

 private:

  static const Int_t fNMaxBinsCentrality = 20;
  static const Int_t fNMaxBinsPt = 10;

  AliAODEvent *fAOD;
  TClonesArray *fTracksCentralBarrel, *fTracksMuonArm;
  AliAODTrack *fTrackCB, *fTrackMA;
  
  Int_t fFilterBitCentralBarrel;
  Double_t fMaxEtaCentralBarrel;

  Double_t fMaxChi2Muon, fMinRAbsMuon, fMaxRAbsMuon;
  Short_t fTriggerMatchLevelMuon;

  TString fTriggerWord;
  Bool_t fIsTriggerSet;

  Int_t fNbinsCent, fNbinsPt;
  
  TAxis *fCentAxis, *fMultAxis, *fPtAxis;

  TH1D *fHistDeltaPhi[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt], *fHistDeltaPhiMix[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt];
  TH2D *fHistNTracksCB_vs_NTracksMA[fNMaxBinsCentrality];

  TH1D *fHistV0Multiplicity, *fHistITSMultiplicity;
  TH1D *fHistCentrality;

  TString fCentMethod;

  Int_t fEvt;

  TList *fOutputList;

  AliAnalysisTaskMuonHadronCorrelations(const AliAnalysisTaskMuonHadronCorrelations&);//not implimented
  AliAnalysisTaskMuonHadronCorrelations& operator=(const AliAnalysisTaskMuonHadronCorrelations&);//not implimnted
  
  ClassDef(AliAnalysisTaskMuonHadronCorrelations, 1)  // example of analysis

};

//====================================================================================================================================================

#endif
