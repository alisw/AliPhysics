#ifndef AliMuonForwardTrackAnalysis_H
#define AliMuonForwardTrackAnalysis_H 

//====================================================================================================================================================
//
//      Class for the analysis of the ALICE muon forward tracks (MUON + MFT)
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TObject.h"
#include "TClonesArray.h"
#include "AliMuonForwardTrack.h"
#include "AliMuonForwardTrackPair.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "AliLog.h"
#include "TFile.h"
#include "TParticle.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "TAxis.h"

//====================================================================================================================================================

class AliMuonForwardTrackAnalysis : public TObject {

public:

  enum {kNoOption, kResonanceOnly, kCharmOnly, kBeautyOnly, kBackground1mu, kBackground2mu, kNoResonances};
  enum {kSingleEvent, kMixedEvent};
  
  AliMuonForwardTrackAnalysis();
  
  virtual ~AliMuonForwardTrackAnalysis() {};  // destructor

  Bool_t Init(Char_t *inputFileName);
  Bool_t LoadNextEvent();
  void Terminate(Char_t *outputFileName);

  void BookHistos();

  void SetInputDir(Char_t *inputDir)   { fInputDir  = inputDir; }
  void SetOutputDir(Char_t *outputDir) { fOutputDir = outputDir; }

  Int_t GetNTracksAnalyzed() { return fNTracksAnalyzed; }

  void SetMassRange(Double_t min, Double_t max) { fMassMin=min; fMassMax=max; }
  void SetTrueMass(Double_t mass) { fTrueMass = mass; }
  void SetSingleMuonAnalysis(Bool_t singleMuonAnalysis) { fSingleMuonAnalysis = singleMuonAnalysis; }
  void SetMuonPairAnalysis(Bool_t muonPairAnalysis) { fMuonPairAnalysis = muonPairAnalysis; }
  void SetMatchTrigger(Int_t triggerLevel) { fTriggerLevel = triggerLevel; }

  Bool_t AnalyzeSingleMuon();
  Bool_t AnalyzeMuonPair(Int_t opt);
  void BuildMuonPairs();
  void BuildMuonPairsMix();

  void EvalDimuonVtxResolution(Bool_t eval) { fEvalDimuonVtxResolution = eval; }

  Bool_t PassedCutSingleMuon(AliMuonForwardTrack *track);
  Bool_t PassedCutMuonPair(AliMuonForwardTrackPair *pair);

  void SetVertResMC(Double_t xRes, Double_t yRes, Double_t zRes) { fXVertResMC=xRes; fYVertResMC=yRes; fZVertResMC=zRes; }

  void SetOption(Int_t option) { fOption = option; }

  void SetMaxNWrongClustersMC(Int_t nClusters) { fMaxNWrongClustersMC = nClusters; }
  void SetMinPtSingleMuons(Double_t ptMin) { fMinPtSingleMuons = ptMin; }
  void SetEtaRangeSingleMuons(Double_t min, Double_t max) { fMinEtaSingleMuons=min; fMaxEtaSingleMuons=max; }
  void SetMaxChi2SingleMuons(Double_t chi2Max) { fMaxChi2SingleMuons = chi2Max; }
  void SetMaxOffsetSingleMuons(Double_t offsetMax) { fMaxOffsetSingleMuons = offsetMax; }
  void CorrelateCutOnOffsetChi2(Bool_t option) { fCorrelateCutOnOffsetChi2 = option; }

  void SetMaxWOffsetMuonPairsAtPrimaryVtx(Double_t wOffsetMax) { fMaxWOffsetMuonPairsAtPrimaryVtx = wOffsetMax; }
  void SetMaxWOffsetMuonPairsAtPCA(Double_t wOffsetMax) { fMaxWOffsetMuonPairsAtPCA = wOffsetMax; }
  void SetMaxDistancePrimaryVtxPCA(Double_t distanceMax) { fMaxDistancePrimaryVtxPCA = distanceMax; }
  void SetMinPCAQuality(Double_t minQuality) { fMinPCAQuality = minQuality; }

  void ReadEvents(Int_t firstEv, Int_t lastEv) { fFirstEvent = firstEv; fLastEvent = lastEv; }
  Int_t GetFirstEvent() { return fFirstEvent; }
  Int_t GetLastEvent()  { return fLastEvent; }

  void UseBransonForCut(Bool_t useBranson) { fUseBransonForCut = useBranson; }
  void UseBransonForKinematics(Bool_t useBranson) { fUseBransonForKinematics = useBranson; }

  Double_t GetPseudoProperDecayLength(AliMuonForwardTrackPair *pair, Double_t trueMass);

  void SetNEventsToMix(Int_t nEvents) { fNEventsToMix = nEvents; if (0<fNEventsToMix && fNEventsToMix<100) fMixing = kTRUE; }

private:

  TString fInputDir, fOutputDir;

  TTree *fInputTreeWithBranson, *fInputTreeWithoutBranson;   //!

  TClonesArray *fMuonForwardTracksWithBranson,    *fMuonForwardTracksWithBransonMix,    *fMuonForwardTrackPairsWithBranson;                //!
  TClonesArray *fMuonForwardTracksWithoutBranson, *fMuonForwardTracksWithoutBransonMix, *fMuonForwardTrackPairsWithoutBranson;             //!
  AliMuonForwardTrack *fMFTTrackWithBranson, *fMFTTrackWithoutBranson, *fMFTTrack;                   //!
  AliMuonForwardTrackPair *fMFTTrackPairWithBranson, *fMFTTrackPairWithoutBranson, *fMFTTrackPair;   //!
  TParticle *fMCRefTrack;                                                                            //!

  Int_t fEv, fEvMix, fFirstEvent, fLastEvent, fNTracksOfEvent, fNTracksAnalyzedOfEvent, fNTracksAnalyzed, fNPairsOfEvent, fNPairsAnalyzedOfEvent;
  Int_t fNTracksAnalyzedOfEventAfterCut, fNPairsAnalyzedOfEventAfterCut;
  
  TH3D *fHistXOffsetSingleMuonsVsEtaVsP,  *fHistYOffsetSingleMuonsVsEtaVsP,  *fHistOffsetSingleMuonsVsEtaVsP,  *fHistWOffsetSingleMuonsVsEtaVsP;      //!
  TH3D *fHistXOffsetSingleMuonsVsEtaVsPt, *fHistYOffsetSingleMuonsVsEtaVsPt, *fHistOffsetSingleMuonsVsEtaVsPt, *fHistWOffsetSingleMuonsVsEtaVsPt;     //!
  TH3D *fHistXErrorSingleMuonsVsEtaVsP,   *fHistYErrorSingleMuonsVsEtaVsP;                                                                  //!
  TH3D *fHistXErrorSingleMuonsVsEtaVsPt,  *fHistYErrorSingleMuonsVsEtaVsPt;                                                                 //!
  TH1D *fHistZOriginSingleMuonsMC;
  
  TH2D *fHistZROriginSingleMuonsMC, *fHistSingleMuonsPtRapidity, *fHistSingleMuonsOffsetChi2;   //!

  TH2D *fHistMassMuonPairsMCVsPt;         //!
  TH2D *fHistDimuonVtxResolutionXVsPt;    //!
  TH2D *fHistDimuonVtxResolutionYVsPt;    //!
  TH2D *fHistDimuonVtxResolutionZVsPt;    //!

  TH2D *fHistRapidityPtMuonPairs[2];              //!
  TH2D *fHistMassMuonPairsVsPt[2];                //!
  TH2D *fHistMassMuonPairsWithoutMFTVsPt[2];      //!
  TH2D *fHistMassMuonPairsVsPtLSp[2];             //!
  TH2D *fHistMassMuonPairsWithoutMFTVsPtLSp[2];   //!
  TH2D *fHistMassMuonPairsVsPtLSm[2];             //!
  TH2D *fHistMassMuonPairsWithoutMFTVsPtLSm[2];   //!

  TH2D *fHistWOffsetMuonPairsAtPrimaryVtxVsPt[2];    //!
  TH2D *fHistWOffsetMuonPairsAtPCAVsPt[2];           //!
  TH2D *fHistDistancePrimaryVtxPCAVsPt[2];           //!
  TH2D *fHistPCAQualityVsPt[2];                      //!
  TH2D *fHistPseudoProperDecayLengthVsPt[2];         //!

  Bool_t fEvalDimuonVtxResolution;

  Double_t fTrueMass, fMassMin, fMassMax;

  Bool_t fSingleMuonAnalysis, fMuonPairAnalysis;
  Int_t fOption, fTriggerLevel;

  Double_t fXVertResMC, fYVertResMC, fZVertResMC;
  Double_t fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ;
  Int_t fMaxNWrongClustersMC;
  Double_t fMinPtSingleMuons, fMinEtaSingleMuons, fMaxEtaSingleMuons;

  Bool_t fUseBransonForCut, fUseBransonForKinematics, fCorrelateCutOnOffsetChi2;

  Double_t fMaxChi2SingleMuons, fMaxOffsetSingleMuons;
  Double_t fMaxWOffsetMuonPairsAtPrimaryVtx, fMaxWOffsetMuonPairsAtPCA, fMaxDistancePrimaryVtxPCA, fMinPCAQuality;

  Bool_t fMixing;
  Int_t fNEventsToMix;

  ClassDef(AliMuonForwardTrackAnalysis, 1)

};

//====================================================================================================================================================
	
#endif

