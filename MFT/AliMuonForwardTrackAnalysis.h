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
#include "TMatrixD.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
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

  enum {kNoOption, kOpenFlavor, kResonanceOnly};
  
  AliMuonForwardTrackAnalysis();
  
  virtual ~AliMuonForwardTrackAnalysis() {};  // destructor

  Bool_t Init(Char_t *inputFileName);
  Bool_t LoadNextEvent();
  void Terminate(Char_t *outputFileName);

  void BookHistos();

  void SetInputDir(Char_t *inputDir)   { fInputDir  = inputDir; }
  void SetOutputDir(Char_t *outputDir) { fOutputDir = outputDir; }

  Int_t GetNTracksAnalyzed() { return fNTracksAnalyzed; }

  void SetMassRange(Int_t nBins, Double_t min, Double_t max) { fNMassBins=nBins; fMassMin=min; fMassMax=max; }
  void SetPtDimuRange(Int_t nBins, Double_t min, Double_t max) { fNPtDimuBins=TMath::Min(nBins,fNMaxPtBinsDimuons) ; fPtDimuMin=min; fPtDimuMax=max; }
  void SetSingleMuonAnalysis(Bool_t singleMuonAnalysis) { fSingleMuonAnalysis = singleMuonAnalysis; }
  void SetMuonPairAnalysis(Bool_t muonPairAnalysis) { fMuonPairAnalysis = muonPairAnalysis; }
  void SetMatchTrigger(Bool_t matchTrigger) { fMatchTrigger = matchTrigger; }

  Bool_t AnalyzeSingleMuon();
  Bool_t AnalyzeMuonPair();
  void BuildMuonPairs();

  Bool_t PassedCutSingleMuon(AliMuonForwardTrack *track);
  Bool_t PassedCutMuonPair(AliMuonForwardTrackPair *pair);

  void SetVertResMC(Double_t xRes, Double_t yRes, Double_t zRes) { fXVertResMC=xRes; fYVertResMC=yRes; fZVertResMC=zRes; }

  void SetOption(Int_t option) { fOption = option; }
  void SetMaxNWrongClustersMC(Int_t nClusters) { fMaxNWrongClustersMC = nClusters; }
  void SetPtMinSingleMuons(Double_t ptMin) { fPtMinSingleMuons = ptMin; }

  void ReadEvents(Int_t firstEv, Int_t lastEv) { fFirstEvent = firstEv; fLastEvent = lastEv; }
  Int_t GetFirstEvent() { return fFirstEvent; }
  Int_t GetLastEvent()  { return fLastEvent; }

  void UseBransonForCut(Bool_t useBranson) { fUseBransonForCut = useBranson; }
  void UseBransonForKinematics(Bool_t useBranson) { fUseBransonForKinematics = useBranson; }
  void UseCutOnOffsetChi2(Bool_t useCut) { fCutOnOffsetChi2 = useCut; }

private:

  static const Int_t fNMaxPtBinsDimuons = 50;

  TString fInputDir, fOutputDir;

  TTree *fInputTreeWithBranson, *fInputTreeWithoutBranson;  //!

  TClonesArray *fMuonForwardTracksWithBranson,    *fMuonForwardTrackPairsWithBranson;                //!
  TClonesArray *fMuonForwardTracksWithoutBranson, *fMuonForwardTrackPairsWithoutBranson;             //!
  AliMuonForwardTrack *fMFTTrackWithBranson, *fMFTTrackWithoutBranson, *fMFTTrack;                   //!
  AliMuonForwardTrackPair *fMFTTrackPairWithBranson, *fMFTTrackPairWithoutBranson, *fMFTTrackPair;   //!
  TParticle *fMCRefTrack;                                                                            //!

  Int_t fEv, fFirstEvent, fLastEvent, fNTracksOfEvent, fNTracksAnalyzedOfEvent, fNTracksAnalyzed, fNPairsOfEvent, fNPairsAnalyzedOfEvent;
  
  TH1D *fHistOffsetSingleMuonsX, *fHistOffsetSingleMuonsY, *fHistOffsetSingleMuons, *fHistWOffsetSingleMuons;      //!
  TH1D *fHistErrorSingleMuonsX, *fHistErrorSingleMuonsY;                                                           //!
  TH2D *fHistSingleMuonsPtRapidity, *fHistSingleMuonsOffsetChi2;                                                   //!
  TGraph *fGraphSingleMuonsOffsetChi2;                                                                             //!

  TH1D *fHistWOffsetMuonPairs[fNMaxPtBinsDimuons+1];          //!
  TH1D *fHistMassMuonPairs[fNMaxPtBinsDimuons+1];             //!
  TH1D *fHistMassMuonPairsWithoutMFT[fNMaxPtBinsDimuons+1];   //!
  TH1D *fHistMassMuonPairsMC[fNMaxPtBinsDimuons+1];           //!
  TH2D *fHistRapidityPtMuonPairsMC;                           //!
 
  Int_t fNMassBins, fNPtDimuBins;
  Double_t fMassMin, fMassMax, fPtDimuMin, fPtDimuMax;
  TAxis *fPtAxisDimuons;

  Bool_t fSingleMuonAnalysis, fMuonPairAnalysis, fMatchTrigger;
  Int_t fOption;

  Double_t fXVertResMC, fYVertResMC, fZVertResMC;
  Double_t fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ;
  Int_t fMaxNWrongClustersMC;
  Double_t fPtMinSingleMuons;

  Bool_t fUseBransonForCut, fUseBransonForKinematics, fCutOnOffsetChi2;

  Double_t fCenterOffset, fCenterChi2, fScaleOffset, fScaleChi2, fRadiusCut;

  ClassDef(AliMuonForwardTrackAnalysis, 1)

};

//====================================================================================================================================================
	
#endif

