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
#include "AliLog.h"
#include "TFile.h"
#include "TParticle.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "TGeoManager.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "AliMuonForwardTrackAnalysis.h"

ClassImp(AliMuonForwardTrackAnalysis)

//====================================================================================================================================================

AliMuonForwardTrackAnalysis::AliMuonForwardTrackAnalysis():
  TObject(),
  fInputDir(0),
  fOutputDir(0),
  fInputTreeWithBranson(0x0),
  fInputTreeWithoutBranson(0x0),
  fMuonForwardTracksWithBranson(0),
  fMuonForwardTrackPairsWithBranson(0),
  fMuonForwardTracksWithoutBranson(0),
  fMuonForwardTrackPairsWithoutBranson(0),
  fMFTTrackWithBranson(0),
  fMFTTrackWithoutBranson(0),
  fMFTTrack(0),
  fMFTTrackPairWithBranson(0),
  fMFTTrackPairWithoutBranson(0),
  fMFTTrackPair(0),
  fMCRefTrack(0),
  fEv(0),
  fFirstEvent(-1),
  fLastEvent(-1),
  fNTracksOfEvent(0),
  fNTracksAnalyzedOfEvent(0),
  fNTracksAnalyzed(0),
  fNPairsOfEvent(0),
  fNPairsAnalyzedOfEvent(0),
  fHistOffsetSingleMuonsX(0x0),
  fHistOffsetSingleMuonsY(0x0),
  fHistOffsetSingleMuons(0x0),
  fHistWOffsetSingleMuons(0x0),
  fHistErrorSingleMuonsX(0x0),
  fHistErrorSingleMuonsY(0x0),
  fHistSingleMuonsPtRapidity(0x0),
  fHistSingleMuonsOffsetChi2(0x0),
  fGraphSingleMuonsOffsetChi2(0x0),
  fHistRapidityPtMuonPairsMC(0x0),
  fNMassBins(10),
  fNPtDimuBins(1000),
  fMassMin(0),
  fMassMax(10),
  fPtDimuMin(0),
  fPtDimuMax(5),
  fPtAxisDimuons(0),
  fSingleMuonAnalysis(1),
  fMuonPairAnalysis(1),
  fMatchTrigger(0),
  fOption(0),
  fXVertResMC(50.e-4),
  fYVertResMC(50.e-4),
  fZVertResMC(50.e-4),
  fPrimaryVtxX(0.),
  fPrimaryVtxY(0.),
  fPrimaryVtxZ(0.),
  fMaxNWrongClustersMC(999),
  fPtMinSingleMuons(0),
  fUseBransonForCut(kFALSE),
  fUseBransonForKinematics(kFALSE),
  fCutOnOffsetChi2(kFALSE),
  fCenterOffset(0.), 
  fCenterChi2(1.), 
  fScaleOffset(250.), 
  fScaleChi2(4.5),
  fRadiusCut(1.)
{

  // default constructor

  for (Int_t iPtBin=0; iPtBin<fNMaxPtBinsDimuons+1; iPtBin++) {
    fHistWOffsetMuonPairs[iPtBin]        = NULL;
    fHistMassMuonPairs[iPtBin]           = NULL;
    fHistMassMuonPairsWithoutMFT[iPtBin] = NULL;
    fHistMassMuonPairsMC[iPtBin]         = NULL;
  }

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::Init(Char_t *inputFileName) {

  fPtAxisDimuons = new TAxis(fNPtDimuBins, fPtDimuMin, fPtDimuMax);

  BookHistos();

  TFile *inputFile = new TFile(Form("%s/%s",fInputDir.Data(),inputFileName));
  if (!inputFile || !inputFile->IsOpen()) {
    AliError(Form("Error opening file %s", inputFileName));
    return kFALSE;
  }
  fInputTreeWithBranson = (TTree*) inputFile->Get("AliMuonForwardTracksWithBranson");
  if (!fInputTreeWithBranson) {
    AliError("Error reading input tree");
    return kFALSE;
  }
  fInputTreeWithoutBranson = (TTree*) inputFile->Get("AliMuonForwardTracksWithoutBranson");
  if (!fInputTreeWithoutBranson) {
    AliError("Error reading input tree");
    return kFALSE;
  }

  if (fFirstEvent<0 || fLastEvent<0 || fFirstEvent>fLastEvent || fFirstEvent>=fInputTreeWithBranson->GetEntries()) {
    fFirstEvent = 0;
    fLastEvent  = fInputTreeWithBranson->GetEntries()-1;
  }
  else {
    fLastEvent = TMath::Min(fLastEvent, Int_t(fInputTreeWithBranson->GetEntries()-1));
  }

  AliInfo(Form("Analysing events %d to %d", fFirstEvent, fLastEvent));

  fMuonForwardTracksWithBranson = new TClonesArray("AliMuonForwardTrack");
  fInputTreeWithBranson->SetBranchAddress("tracks", &fMuonForwardTracksWithBranson);  

  fMuonForwardTracksWithoutBranson = new TClonesArray("AliMuonForwardTrack");
  fInputTreeWithoutBranson->SetBranchAddress("tracks", &fMuonForwardTracksWithoutBranson);  

  TGeoManager::Import(Form("%s/geometry.root",fInputDir.Data()));

  AliMUONTrackExtrap::SetField();

  fMuonForwardTrackPairsWithBranson    = new TClonesArray("AliMuonForwardTrackPair");
  fMuonForwardTrackPairsWithoutBranson = new TClonesArray("AliMuonForwardTrackPair");

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::LoadNextEvent() {

  if (fEv>fLastEvent) return kFALSE;
  if (fEv<fFirstEvent) { fEv++; return kTRUE; }
  fMuonForwardTracksWithBranson -> Delete();
  fMuonForwardTracksWithoutBranson -> Delete();
  fInputTreeWithBranson->GetEvent(fEv);
  fInputTreeWithoutBranson->GetEvent(fEv);
  AliInfo(Form("**** analyzing event # %4d (%3d tracks) ****", fEv, fMuonForwardTracksWithBranson->GetEntries()));

  fPrimaryVtxX = gRandom->Gaus(0., fXVertResMC);
  fPrimaryVtxY = gRandom->Gaus(0., fYVertResMC);
  fPrimaryVtxZ = gRandom->Gaus(0., fZVertResMC);

  if (fSingleMuonAnalysis) {
    fNTracksAnalyzedOfEvent = 0;
    fNTracksOfEvent = fMuonForwardTracksWithBranson->GetEntries();
    while (AnalyzeSingleMuon()) continue;
  }
  
  if (fMuonPairAnalysis) {
    if (fMuonForwardTrackPairsWithBranson) {
      fMuonForwardTrackPairsWithBranson->Delete();
      fMuonForwardTrackPairsWithoutBranson->Delete();
    }
    BuildMuonPairs();
    fNPairsAnalyzedOfEvent = 0;
    fNPairsOfEvent = fMuonForwardTrackPairsWithBranson->GetEntries();
    while (AnalyzeMuonPair()) continue;
  }

  fEv++;
  
  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeSingleMuon() {

  if (fNTracksAnalyzedOfEvent>=fNTracksOfEvent) return kFALSE;

  fMFTTrackWithBranson    = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson->At(fNTracksAnalyzedOfEvent);
  fMFTTrackWithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson->At(fNTracksAnalyzedOfEvent);

  fNTracksAnalyzedOfEvent++;

  Bool_t passedCut = kFALSE;
  if (fUseBransonForCut) passedCut = PassedCutSingleMuon(fMFTTrackWithBranson);
  else passedCut = PassedCutSingleMuon(fMFTTrackWithoutBranson);
  if (!passedCut) return kTRUE;

  if (fUseBransonForKinematics) fMFTTrack = fMFTTrackWithBranson;
  else fMFTTrack = fMFTTrackWithoutBranson;
  if (fMatchTrigger && !fMFTTrack->GetMatchTrigger()) return kTRUE;
  fMCRefTrack = fMFTTrack->GetMCTrackRef();

  if (!fMCRefTrack) return kTRUE;
  if (fMFTTrack->GetNWrongClustersMC()>fMaxNWrongClustersMC) return kTRUE;

  AliMUONTrackParam *param = fMFTTrack->GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, fPrimaryVtxZ);

  TLorentzVector pMu;
  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(param->P()*param->P() + mMu*mMu);
  pMu.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  TMatrixD cov(5,5);
  cov = param->GetCovariances();

  fHistErrorSingleMuonsX -> Fill(1.e4*TMath::Sqrt(cov(0,0)));
  fHistErrorSingleMuonsY -> Fill(1.e4*TMath::Sqrt(cov(2,2)));

  Double_t dX = fMFTTrack->GetOffsetX(fPrimaryVtxX, fPrimaryVtxZ);
  Double_t dY = fMFTTrack->GetOffsetY(fPrimaryVtxY, fPrimaryVtxZ);
  
  Double_t offset = fMFTTrack->GetOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);
  Double_t weightedOffset = fMFTTrack->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);

  //  AliDebug(2, Form("pdg code = %d\n", fMCRefTrack->GetPdgCode()));

  fHistOffsetSingleMuonsX -> Fill(1.e4*dX);
  fHistOffsetSingleMuonsY -> Fill(1.e4*dY);

  fHistSingleMuonsPtRapidity  -> Fill(pMu.Rapidity(), pMu.Pt());
  fHistOffsetSingleMuons      -> Fill(1.e4*offset);
  fHistWOffsetSingleMuons     -> Fill(weightedOffset);
  Double_t chi2OverNdf = fMFTTrack->GetGlobalChi2()/Double_t(fMFTTrack->GetNMFTClusters()+fMFTTrack->GetNMUONClusters());
  fHistSingleMuonsOffsetChi2  -> Fill(1.e4*offset, chi2OverNdf);
  fGraphSingleMuonsOffsetChi2 -> SetPoint(fGraphSingleMuonsOffsetChi2->GetN(),1.e4*offset, chi2OverNdf);

  fNTracksAnalyzed++;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeMuonPair() {

  if (fNPairsAnalyzedOfEvent>=fNPairsOfEvent) return kFALSE;

  fMFTTrackPairWithBranson    = (AliMuonForwardTrackPair*) fMuonForwardTrackPairsWithBranson->At(fNPairsAnalyzedOfEvent);
  fMFTTrackPairWithoutBranson = (AliMuonForwardTrackPair*) fMuonForwardTrackPairsWithoutBranson->At(fNPairsAnalyzedOfEvent);

  fNPairsAnalyzedOfEvent++;

  Bool_t passedCut = kFALSE;
  if (fUseBransonForCut) passedCut = PassedCutMuonPair(fMFTTrackPairWithBranson);
  else passedCut = PassedCutMuonPair(fMFTTrackPairWithoutBranson);

  if (!passedCut) return kTRUE;

  if (fUseBransonForKinematics) fMFTTrackPair = fMFTTrackPairWithBranson;
  else fMFTTrackPair = fMFTTrackPairWithoutBranson;

  if ( fMFTTrackPair->GetTrack(0)->GetNWrongClustersMC()>fMaxNWrongClustersMC || 
       fMFTTrackPair->GetTrack(1)->GetNWrongClustersMC()>fMaxNWrongClustersMC ) return kTRUE;

  if (fOption==kResonanceOnly && !fMFTTrackPair->IsResonance()) return kTRUE;

  Int_t ptBin = fPtAxisDimuons->FindBin(fMFTTrackPair->GetPtMC());

  if (1<=ptBin && ptBin<=fNPtDimuBins) {
    fHistMassMuonPairs[ptBin]           -> Fill(fMFTTrackPair->GetMass(fPrimaryVtxZ));
    fHistWOffsetMuonPairs[ptBin]        -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
    fHistMassMuonPairsWithoutMFT[ptBin] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
    fHistMassMuonPairsMC[ptBin]         -> Fill(fMFTTrackPair->GetMassMC());
  }
  fHistMassMuonPairs[0]           -> Fill(fMFTTrackPair->GetMass(fPrimaryVtxZ));
  fHistWOffsetMuonPairs[0]        -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
  fHistMassMuonPairsWithoutMFT[0] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
  fHistMassMuonPairsMC[0]         -> Fill(fMFTTrackPair->GetMassMC());

  fHistRapidityPtMuonPairsMC          -> Fill(fMFTTrackPair->GetRapidityMC(), fMFTTrackPair->GetPtMC());

  AliDebug(1, Form("mass = %f   MC = %f", fMFTTrackPair->GetMass(fPrimaryVtxZ), fMFTTrackPair->GetMassMC()));

  return kTRUE;

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BuildMuonPairs() {

  for (Int_t iTrack=0; iTrack<fMuonForwardTracksWithBranson->GetEntries(); iTrack++) {
    for (Int_t jTrack=0; jTrack<iTrack; jTrack++) {
    
      AliMuonForwardTrack *track0_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson->At(iTrack);
      AliMuonForwardTrack *track1_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson->At(jTrack);

      AliMuonForwardTrack *track0_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson->At(iTrack);
      AliMuonForwardTrack *track1_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson->At(jTrack);

      if (fMatchTrigger) if (!track0_WithBranson->GetMatchTrigger() || !track1_WithBranson->GetMatchTrigger()) continue;
      if (!track0_WithBranson->GetMCTrackRef() || !track1_WithBranson->GetMCTrackRef()) continue;

      AliMuonForwardTrackPair *trackPairWithBranson    = new AliMuonForwardTrackPair(track0_WithBranson, track1_WithBranson);
      AliMuonForwardTrackPair *trackPairWithoutBranson = new AliMuonForwardTrackPair(track0_WithoutBranson, track1_WithoutBranson);
      if (fOption==kResonanceOnly && !trackPairWithBranson->IsResonance()) {
	delete trackPairWithBranson;
	delete trackPairWithoutBranson;
	continue;
      }
      new ((*fMuonForwardTrackPairsWithBranson)[fMuonForwardTrackPairsWithBranson->GetEntries()]) AliMuonForwardTrackPair(*trackPairWithBranson);
      new ((*fMuonForwardTrackPairsWithoutBranson)[fMuonForwardTrackPairsWithoutBranson->GetEntries()]) AliMuonForwardTrackPair(*trackPairWithoutBranson);
    }
  }

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::PassedCutSingleMuon(AliMuonForwardTrack *track) {

  AliMUONTrackParam *param = track->GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, fPrimaryVtxZ);

  if (track->Pt()<fPtMinSingleMuons) return kFALSE;
  
  if (fCutOnOffsetChi2) {
    Double_t offset = 1.e4*track->GetOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);
    Double_t chi2OverNdf = track->GetGlobalChi2() / Double_t(track->GetNMFTClusters()+track->GetNMUONClusters()); 
    offset /= fScaleOffset;
    chi2OverNdf /= fScaleChi2;
    offset -= fCenterOffset;
    chi2OverNdf -= fCenterChi2;
    //    printf("cut on offset and chi2: returning %d\n", TMath::Sqrt(offset*offset + chi2OverNdf*chi2OverNdf)>fRadiusCut);
    if (TMath::Sqrt(offset*offset + chi2OverNdf*chi2OverNdf) > fRadiusCut) return kFALSE;
  }

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::PassedCutMuonPair(AliMuonForwardTrackPair *pair) {

  return PassedCutSingleMuon(pair->GetTrack(0)) && PassedCutSingleMuon(pair->GetTrack(1));

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::Terminate(Char_t *outputFileName) {

  TFile *fileOut = new TFile(Form("%s/%s",fOutputDir.Data(),outputFileName), "recreate");

  printf("Writing output objects to file %s\n", fileOut->GetName());

  fHistOffsetSingleMuonsX -> Write();
  fHistOffsetSingleMuonsY -> Write();
  fHistOffsetSingleMuons  -> Write();
  fHistWOffsetSingleMuons -> Write();
  fHistErrorSingleMuonsX  -> Write();
  fHistErrorSingleMuonsY  -> Write();

  fHistSingleMuonsPtRapidity -> Write();
  fHistSingleMuonsOffsetChi2 -> Write();

  fGraphSingleMuonsOffsetChi2 -> Write();

  for (Int_t iPtBin=0; iPtBin<fNPtDimuBins+1; iPtBin++) {
    fHistWOffsetMuonPairs[iPtBin]        -> Write();
    fHistMassMuonPairs[iPtBin]           -> Write();
    fHistMassMuonPairsWithoutMFT[iPtBin] -> Write();
    fHistMassMuonPairsMC[iPtBin]         -> Write();
  }

  fHistRapidityPtMuonPairsMC -> Write();

  fileOut -> Close();

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BookHistos() {

  fHistOffsetSingleMuonsX = new TH1D("fHistOffsetSingleMuonsX", "Offset for single muons along X",  200, -1000, 1000);
  fHistOffsetSingleMuonsY = new TH1D("fHistOffsetSingleMuonsY", "Offset for single muons along Y",  200, -1000, 1000);
  fHistErrorSingleMuonsX  = new TH1D("fHistErrorSingleMuonsX",  "Coordinate Error for single muons along X",  200, 0, 1000);
  fHistErrorSingleMuonsY  = new TH1D("fHistErrorSingleMuonsY",  "Coordinate Error for single muons along Y",  200, 0, 1000);
  fHistOffsetSingleMuons  = new TH1D("fHistOffsetSingleMuons",  "Offset for single muons",          200, 0, 2000);
  fHistWOffsetSingleMuons = new TH1D("fHistWOffsetSingleMuons", "Weighted Offset for single muons", 300, 0, 15);  

  fHistSingleMuonsPtRapidity = new TH2D("fHistSingleMuonsPtRapidity", "Phase Space for single muons", 10, -4, -2.5, 10, 0.5, 5.5);
  fHistSingleMuonsOffsetChi2 = new TH2D("fHistSingleMuonsOffsetChi2", "Offset vs #chi^{2}/ndf for single muons", 400, 0, 4000, 100, 0, 20);

  fHistOffsetSingleMuonsX -> SetXTitle("Offset(X)  [#mum]");
  fHistOffsetSingleMuonsY -> SetXTitle("Offset(Y)  [#mum]");
  fHistErrorSingleMuonsX  -> SetXTitle("Err. on track position at z_{vtx} (X)  [#mum]");
  fHistErrorSingleMuonsY  -> SetXTitle("Err. on track position at z_{vtx} (Y)  [#mum]");
  fHistOffsetSingleMuons  -> SetXTitle("Offset  [#mum]");
  fHistWOffsetSingleMuons -> SetXTitle("Weighted Offset");

  fHistSingleMuonsPtRapidity -> SetXTitle("y^{#mu}");
  fHistSingleMuonsPtRapidity -> SetYTitle("p_{T}^{#mu}  [GeV/c]");
  fHistSingleMuonsOffsetChi2 -> SetXTitle("Offset  [#mum]");
  fHistSingleMuonsOffsetChi2 -> SetYTitle("#chi^{2}/ndf");

  fHistOffsetSingleMuonsX -> Sumw2();
  fHistOffsetSingleMuonsY -> Sumw2();
  fHistErrorSingleMuonsX  -> Sumw2();
  fHistErrorSingleMuonsY  -> Sumw2();
  fHistOffsetSingleMuons  -> Sumw2();
  fHistWOffsetSingleMuons -> Sumw2();

  fHistSingleMuonsPtRapidity -> Sumw2();
  fHistSingleMuonsOffsetChi2 -> Sumw2();
    
  //--------------------------------------------

  fGraphSingleMuonsOffsetChi2 = new TGraph();
  fGraphSingleMuonsOffsetChi2 -> SetName("fGraphSingleMuonsOffsetChi2");
  fGraphSingleMuonsOffsetChi2 -> SetTitle("fGraphSingleMuonsOffsetChi2");

  //--------------------------------------------

  for (Int_t iPtBin=0; iPtBin<=fNPtDimuBins+1; iPtBin++) {

    if (!iPtBin) {
      fHistWOffsetMuonPairs[iPtBin]        = new TH1D(Form("fHistWOffsetMuonPairs_%d",iPtBin),        
						      "Weighted Offset for Muon Pairs (All p_{T}^{#mu#mu})",
						      300, 0, 60);
      fHistMassMuonPairs[iPtBin]	   = new TH1D(Form("fHistMassMuonPairs_%d",iPtBin),           
						      "Dimuon Mass (MUON+MFT) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsWithoutMFT[iPtBin] = new TH1D(Form("fHistMassMuonPairsWithoutMFT_%d",iPtBin), 
						      "Dimuon Mass (MUON only) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsMC[iPtBin]         = new TH1D(Form("fHistMassMuonPairsMC_%d",iPtBin),         
						      "Dimuon Mass (MC) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
    }

    else {
      Double_t ptMin = fPtAxisDimuons->GetBinLowEdge(iPtBin);
      Double_t ptMax = fPtAxisDimuons->GetBinUpEdge(iPtBin);
      fHistWOffsetMuonPairs[iPtBin]        = new TH1D(Form("fHistWOffsetMuonPairs_%d",iPtBin),        
						      Form("Weighted Offset for Muon Pairs ( %3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      300, 0, 60);
      fHistMassMuonPairs[iPtBin]	   = new TH1D(Form("fHistMassMuonPairs_%d",iPtBin),           
						      Form("Dimuon Mass (MUON+MFT) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsWithoutMFT[iPtBin] = new TH1D(Form("fHistMassMuonPairsWithoutMFT_%d",iPtBin), 
						      Form("Dimuon Mass (MUON only) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsMC[iPtBin]         = new TH1D(Form("fHistMassMuonPairsMC_%d",iPtBin),         
						      Form("Dimuon Mass (MC) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);      
    }
    
    fHistWOffsetMuonPairs[iPtBin]        -> SetXTitle("Weighted Offset");
    fHistMassMuonPairs[iPtBin]           -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsWithoutMFT[iPtBin] -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsMC[iPtBin]         -> SetXTitle("Mass  [GeV/c^{2}]");    

    fHistWOffsetMuonPairs[iPtBin]        -> Sumw2();
    fHistMassMuonPairs[iPtBin]           -> Sumw2();
    fHistMassMuonPairsWithoutMFT[iPtBin] -> Sumw2();
    fHistMassMuonPairsMC[iPtBin]         -> Sumw2();    

  }
  
  fHistRapidityPtMuonPairsMC   = new TH2D("fHistRapidityPtMuonPairsMC", "Dimuon Phase Space (MC)", 100, -4.5, -2., 100, 0., 10.); 
  fHistRapidityPtMuonPairsMC   -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairsMC   -> SetYTitle("p_{T}  [GeV/c]");
  fHistRapidityPtMuonPairsMC   -> Sumw2();

}

//====================================================================================================================================================

