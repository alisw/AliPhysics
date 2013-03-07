#ifndef AliAnalysisTaskDiMuonCorrelations_H
#define AliAnalysisTaskDiMuonCorrelations_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TH2D.h"

class AliEventPoolManager;

//====================================================================================================================================================

class  AliAnalysisTaskDiMuonCorrelations : public AliAnalysisTaskSE {

 public:
 
  enum {kSingleEvent, kMixedEvent};

  AliAnalysisTaskDiMuonCorrelations();
  AliAnalysisTaskDiMuonCorrelations(const Char_t *name);
  virtual ~AliAnalysisTaskDiMuonCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // ------------- Cuts -----------------

  void SetTriggerMatchLevelMuon(Short_t level) { fTriggerMatchLevelMuon = level; }
  //  void SetMaxChi2Muon(Double_t chi2Max) { fMaxChi2Muon = chi2Max; }
  void SetEtaRangeMuon (Double_t etaMin, Double_t etaMax) { fMinEtaMuon = etaMin; fMaxEtaMuon = etaMax; }

  // ------------- Analysis -------------

  Float_t GetV0Multiplicity();
  Double_t GetITSMultiplicity();
  Bool_t IsTriggerFired();
  TObjArray* GetAcceptedTracksMuonArm(AliAODEvent *aodEvent);
  void SetPtBinning(Int_t nBins, Double_t *limits);
  void SetCentBinning(Int_t nBins, Double_t *limits);
  void SetEtaBinning(Int_t nBins, Double_t *limits);
  void SetCentMethod(const Char_t *method) { fCentMethod = method; }
  void FillHistograms(Int_t centrality, Int_t option);
  Int_t GetCentBin();

 private:

  static const Int_t fNMaxBinsCentrality = 20;
  static const Int_t fNMaxBinsPt = 10;
  static const Int_t fNMaxBinsEta = 10;

  AliAODEvent *fAOD; //!
  AliEventPoolManager *fPoolMgr; //! event pool manager
  AliAODTrack *fMuonTrack[2]; //!
  
  Double_t fMaxChi2Muon, fMinEtaMuon, fMaxEtaMuon;
  Short_t fTriggerMatchLevelMuon;

  Int_t fNbinsCent, fNbinsPt;
  
  TAxis *fCentAxis;
  TAxis *fPtAxis;
  TAxis *fEtaAxis;

  TH1D *fHistDeltaPhi[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
  TH1D *fHistDeltaPhiMix[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
  TH2D *fHistEtaDeltaPhi[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
  TH2D *fHistEtaDeltaPhiMix[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsPt]; //!
  TH2D *fHistNMuons_vs_NMuons[fNMaxBinsCentrality]; //!
  TH2D *fHistNMuons_vs_NMuons_Mixed[fNMaxBinsCentrality]; //!
  TH2D *fHistTracksEtavsEta[fNMaxBinsCentrality]; //!
  TH2D *fHistTracksEtavsEta_Mixed[fNMaxBinsCentrality]; //!
  TH1D *fHistSingleMuonsPt[fNMaxBinsCentrality]; //!
  TH1D *fHistSingleMuonsPt_Mixed[fNMaxBinsCentrality]; //!
  TH2D *fHistSingleMuonsEtaPt[fNMaxBinsCentrality]; //!
  TH2D *fHistSingleMuonsEtaPt_Mixed[fNMaxBinsCentrality]; //!

  TH1D *fHistV0Multiplicity; //!
  TH1D *fHistITSMultiplicity; //!
  TH1D *fHistCentrality; //!
  TH1D *fHistEvStat; //!

  TString fCentMethod;

  TList *fOutputList; //!

  AliAnalysisTaskDiMuonCorrelations(const AliAnalysisTaskDiMuonCorrelations&);//not implimented
  AliAnalysisTaskDiMuonCorrelations& operator=(const AliAnalysisTaskDiMuonCorrelations&);//not implimnted
  
  ClassDef(AliAnalysisTaskDiMuonCorrelations, 2)  // example of analysis

};

//====================================================================================================================================================

#endif
