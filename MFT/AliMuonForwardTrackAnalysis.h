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
  void SetSingleMuonAnalysis(Bool_t singleMuonAnalysis) { fSingleMuonAnalysis = singleMuonAnalysis; }
  void SetMuonPairAnalysis(Bool_t muonPairAnalysis) { fMuonPairAnalysis = muonPairAnalysis; }
  void SetMatchTrigger(Bool_t matchTrigger) { fMatchTrigger = matchTrigger; }

  Bool_t AnalyzeSingleMuon();
  Bool_t AnalyzeMuonPair();
  void BuildMuonPairs();

  void SetVertResMC(Double_t xRes, Double_t yRes, Double_t zRes) { fXVertResMC=xRes; fYVertResMC=yRes; fZVertResMC=zRes; }

  void SetOption(Int_t option) { fOption = option; }
  void SetMaxNWrongClustersMC(Int_t nClusters) { fMaxNWrongClustersMC = nClusters; }
  void SetPtMinSingleMuons(Double_t ptMin) { fPtMinSingleMuons = ptMin; }

  void ReadEvents(Int_t firstEv, Int_t lastEv) { fFirstEvent = firstEv; fLastEvent = lastEv; }
  Int_t GetFirstEvent() { return fFirstEvent; }
  Int_t GetLastEvent()  { return fLastEvent; }

private:

  static const Int_t fNPtBinsOffsetSingleMuons  = 10;
  static const Int_t fNRapBinsOffsetSingleMuons = 10;

  TString fInputDir, fOutputDir;

  TTree *fInputTree;  //!

  TClonesArray *fMuonForwardTracks, *fMuonForwardTrackPairs;  //!
  AliMuonForwardTrack *fMFTTrack;                             //!
  AliMuonForwardTrackPair *fMFTTrackPair;                     //!
  TParticle *fMCRefTrack;                                     //!

  Int_t fEv, fFirstEvent, fLastEvent, fNTracksOfEvent, fNTracksAnalyzedOfEvent, fNTracksAnalyzed, fNPairsOfEvent, fNPairsAnalyzedOfEvent;
  
  TH1D *fHistOffsetSingleMuonsX, *fHistOffsetSingleMuonsY, *fHistOffsetSingleMuons, *fHistWOffsetSingleMuons;      //!
  TH1D *fHistErrorSingleMuonsX, *fHistErrorSingleMuonsY;                                                           //!
  TH2D *fHistOffsetSingleMuonsX_vsPtRapidity, *fHistOffsetSingleMuonsY_vsPtRapidity, *fHistSingleMuonsPtRapidity;  //!
  TH1D *fHistOffsetSingleMuonsX_tmp[fNRapBinsOffsetSingleMuons][fNPtBinsOffsetSingleMuons];			   //!
  TH1D *fHistOffsetSingleMuonsY_tmp[fNRapBinsOffsetSingleMuons][fNPtBinsOffsetSingleMuons];			   //!
  TH1D *fHistWOffsetMuonPairs, *fHistMassMuonPairs, *fHistMassMuonPairsWithoutMFT, *fHistMassMuonPairsMC;          //!
  TH2D *fHistRapidityPtMuonPairsMC;
 
  Int_t fNMassBins;
  Double_t fMassMin, fMassMax;

  Bool_t fSingleMuonAnalysis, fMuonPairAnalysis, fMatchTrigger;
  Int_t fOption;

  Double_t fXVertResMC, fYVertResMC, fZVertResMC;
  Int_t fMaxNWrongClustersMC;
  Double_t fPtMinSingleMuons;

  ClassDef(AliMuonForwardTrackAnalysis, 1)

};

//====================================================================================================================================================
	
#endif

