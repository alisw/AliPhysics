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
#include "AliMuonForwardTrackAnalysis.h"

ClassImp(AliMuonForwardTrackAnalysis)

//====================================================================================================================================================

AliMuonForwardTrackAnalysis::AliMuonForwardTrackAnalysis():
  TObject(),
  fInputDir(0),
  fOutputDir(0),
  fInputTree(0),
  fMuonForwardTracks(0),
  fMuonForwardTrackPairs(0),
  fMFTTrack(0),
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
  fHistOffsetSingleMuonsX_vsPtRapidity(0x0),
  fHistOffsetSingleMuonsY_vsPtRapidity(0x0),
  fHistSingleMuonsPtRapidity(0x0),
  fHistWOffsetMuonPairs(0x0),
  fHistMassMuonPairs(0x0),
  fHistMassMuonPairsWithoutMFT(0x0),
  fHistMassMuonPairsMC(0x0),
  fHistRapidityPtMuonPairsMC(0x0),
  fNMassBins(1000),
  fMassMin(0),
  fMassMax(10),
  fSingleMuonAnalysis(1),
  fMuonPairAnalysis(1),
  fMatchTrigger(0),
  fOption(0),
  fXVertResMC(150.e-4),
  fYVertResMC(150.e-4),
  fZVertResMC(100.e-4),
  fMaxNWrongClustersMC(999),
  fPtMinSingleMuons(0)
{

  // default constructor

  for (Int_t rapBin=0; rapBin<fNRapBinsOffsetSingleMuons; rapBin++) {
    for (Int_t ptBin=0; ptBin<fNPtBinsOffsetSingleMuons; ptBin++) {
      fHistOffsetSingleMuonsX_tmp[rapBin][ptBin] = 0x0;
      fHistOffsetSingleMuonsY_tmp[rapBin][ptBin] = 0x0;
    }
  }

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::Init(Char_t *inputFileName) {

  BookHistos();

  TFile *inputFile = new TFile(Form("%s/%s",fInputDir.Data(),inputFileName));
  if (!inputFile || !inputFile->IsOpen()) {
    AliError(Form("Error opening file %s", inputFileName));
    return kFALSE;
  }
  fInputTree = (TTree*) inputFile->Get("AliMuonForwardTracks");
  if (!fInputTree) {
    AliError("Error reading input tree");
    return kFALSE;
  }

  if (fFirstEvent<0 || fLastEvent<0 || fFirstEvent>fLastEvent || fFirstEvent>=fInputTree->GetEntries()) {
    fFirstEvent = 0;
    fLastEvent  = fInputTree->GetEntries()-1;
  }
  else {
    fLastEvent = TMath::Min(fLastEvent, Int_t(fInputTree->GetEntries()-1));
  }

  AliInfo(Form("Analysing events %d to %d", fFirstEvent, fLastEvent));

  fMuonForwardTracks = new TClonesArray("AliMuonForwardTrack");
  fInputTree->SetBranchAddress("tracks", &fMuonForwardTracks);  

  TGeoManager::Import(Form("%s/geometry.root",fInputDir.Data()));

  AliMUONTrackExtrap::SetField();

  fMuonForwardTrackPairs = new TClonesArray("AliMuonForwardTrackPair");

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::LoadNextEvent() {

  if (fEv>fLastEvent) return kFALSE;
  if (fEv<fFirstEvent) { fEv++; return kTRUE; }
  fMuonForwardTracks -> Clear();
  fInputTree->GetEvent(fEv);
  AliInfo(Form("**** analyzing event # %4d (%3d tracks) ****", fEv, fMuonForwardTracks->GetEntries()));

  if (fSingleMuonAnalysis) {
    fNTracksAnalyzedOfEvent = 0;
    fNTracksOfEvent = fMuonForwardTracks->GetEntries();
    while (AnalyzeSingleMuon()) continue;
  }
  
  if (fMuonPairAnalysis) {
    if (fMuonForwardTrackPairs) {
      fMuonForwardTrackPairs->Clear();
    }
    BuildMuonPairs();
    fNPairsAnalyzedOfEvent = 0;
    fNPairsOfEvent = fMuonForwardTrackPairs->GetEntries();
    while (AnalyzeMuonPair()) continue;
  }

  fEv++;
  
  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeSingleMuon() {

  if (fNTracksAnalyzedOfEvent>=fNTracksOfEvent) return kFALSE;

  fMFTTrack   = (AliMuonForwardTrack*) fMuonForwardTracks->At(fNTracksAnalyzedOfEvent);
  fNTracksAnalyzedOfEvent++;
  if (fMatchTrigger && !fMFTTrack->GetMatchTrigger()) return kTRUE;
  fMCRefTrack = fMFTTrack->GetMCTrackRef();

  if (!fMCRefTrack) return kTRUE;
  if (fMFTTrack->GetNWrongClustersMC()>fMaxNWrongClustersMC) return kTRUE;

  Double_t xOrig=gRandom->Gaus(0., fXVertResMC);
  Double_t yOrig=gRandom->Gaus(0., fXVertResMC);
  Double_t zOrig=gRandom->Gaus(0., fZVertResMC);

  AliMUONTrackParam *param = fMFTTrack->GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, zOrig);

  TLorentzVector pMu;
  Double_t mMu = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t energy = TMath::Sqrt(param->P()*param->P() + mMu*mMu);
  pMu.SetPxPyPzE(param->Px(), param->Py(), param->Pz(), energy);

  if (fMFTTrack->Pt()<fPtMinSingleMuons) return kTRUE;

  TMatrixD cov(5,5);
  cov = param->GetCovariances();

  fHistErrorSingleMuonsX -> Fill(1.e4*TMath::Sqrt(cov(0,0)));
  fHistErrorSingleMuonsY -> Fill(1.e4*TMath::Sqrt(cov(2,2)));

  Double_t dX = fMFTTrack->GetOffsetX(xOrig, zOrig);
  Double_t dY = fMFTTrack->GetOffsetY(yOrig, zOrig);
  
  Double_t offset = fMFTTrack->GetOffset(xOrig, yOrig, zOrig);
  Double_t weightedOffset = fMFTTrack->GetWeightedOffset(xOrig, yOrig, zOrig);

  //  AliDebug(2, Form("pdg code = %d\n", fMCRefTrack->GetPdgCode()));

  fHistOffsetSingleMuonsX -> Fill(1.e4*dX);
  fHistOffsetSingleMuonsY -> Fill(1.e4*dY);
  Int_t rapBin = fHistOffsetSingleMuonsX_vsPtRapidity->GetXaxis()->FindBin(pMu.Rapidity());
  Int_t ptBin  = fHistOffsetSingleMuonsX_vsPtRapidity->GetYaxis()->FindBin(pMu.Pt());
  //  AliDebug(1, Form("pt = %f (%d), rap = %f (%d)\n", pMu.Pt(), pMu.Rapidity(), ptBin, rapBin));
  if (0<rapBin && rapBin<=fNRapBinsOffsetSingleMuons && 0<ptBin && ptBin<=fNPtBinsOffsetSingleMuons) {
    fHistOffsetSingleMuonsX_tmp[rapBin-1][ptBin-1]->Fill(1.e4*dX);
    fHistOffsetSingleMuonsY_tmp[rapBin-1][ptBin-1]->Fill(1.e4*dY);
  }
  fHistSingleMuonsPtRapidity -> Fill(pMu.Rapidity(), pMu.Pt());
  fHistOffsetSingleMuons     -> Fill(1.e4*offset);
  fHistWOffsetSingleMuons    -> Fill(weightedOffset);

  fNTracksAnalyzed++;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeMuonPair() {

  if (fNPairsAnalyzedOfEvent>=fNPairsOfEvent) return kFALSE;

  fMFTTrackPair = (AliMuonForwardTrackPair*) fMuonForwardTrackPairs->At(fNPairsAnalyzedOfEvent);

  if (fOption==kResonanceOnly && !fMFTTrackPair->IsResonance()) {
    fNPairsAnalyzedOfEvent++;
    return kTRUE;
  }

  Double_t xOrig=gRandom->Gaus(0., fXVertResMC);
  Double_t yOrig=gRandom->Gaus(0., fXVertResMC);
  Double_t zOrig=gRandom->Gaus(0., fZVertResMC);
  AliDebug(1, Form("origin = (%f, %f, %f)", xOrig, yOrig, zOrig));

  fHistMassMuonPairs           -> Fill(fMFTTrackPair->GetMass(zOrig));
  fHistWOffsetMuonPairs        -> Fill(fMFTTrackPair->GetWeightedOffset(xOrig, yOrig, zOrig));
  fHistMassMuonPairsWithoutMFT -> Fill(fMFTTrackPair->GetMassWithoutMFT(xOrig, yOrig, zOrig));
  fHistMassMuonPairsMC         -> Fill(fMFTTrackPair->GetMassMC());
  fHistRapidityPtMuonPairsMC   -> Fill(fMFTTrackPair->GetRapidityMC(), fMFTTrackPair->GetPtMC());

  AliDebug(1, Form("mass = %f   MC = %f", fMFTTrackPair->GetMass(zOrig), fMFTTrackPair->GetMassMC()));

  fNPairsAnalyzedOfEvent++;

  return kTRUE;

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BuildMuonPairs() {

  for (Int_t iTrack=0; iTrack<fMuonForwardTracks->GetEntries(); iTrack++) {
    for (Int_t jTrack=0; jTrack<iTrack; jTrack++) {
    
      AliMuonForwardTrack *track0 = (AliMuonForwardTrack*) fMuonForwardTracks->At(iTrack);
      AliMuonForwardTrack *track1 = (AliMuonForwardTrack*) fMuonForwardTracks->At(jTrack);

      if (fMatchTrigger) if (!track0->GetMatchTrigger() || !track1->GetMatchTrigger()) continue;
      if (!track0->GetMCTrackRef() || !track1->GetMCTrackRef()) continue;
      if (track0->GetNWrongClustersMC()>fMaxNWrongClustersMC || track1->GetNWrongClustersMC()>fMaxNWrongClustersMC) continue;
      if (track0->Pt()<fPtMinSingleMuons || track1->Pt()<fPtMinSingleMuons) continue;

      AliMuonForwardTrackPair *trackPair = new AliMuonForwardTrackPair(track0, track1);
      if (fOption==kResonanceOnly && !trackPair->IsResonance()) {
	delete trackPair;
	continue;
      }
      new ((*fMuonForwardTrackPairs)[fMuonForwardTrackPairs->GetEntries()]) AliMuonForwardTrackPair(*trackPair);

    }
  }

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::Terminate(Char_t *outputFileName) {

  for (Int_t rapBin=0; rapBin<fNRapBinsOffsetSingleMuons; rapBin++) {
    for (Int_t ptBin=0; ptBin<fNPtBinsOffsetSingleMuons; ptBin++) {
      Int_t binMin_x = fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->FindBin(-3*fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->GetRMS());
      Int_t binMax_x = fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->FindBin(+3*fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->GetRMS());
      for (Int_t bin=1; bin<=fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->GetNbinsX(); bin++) {
	if (bin<binMin_x || bin>binMax_x) fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->SetBinContent(bin, 0);
      }
      Int_t binMin_y = fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->FindBin(-3*fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->GetRMS());
      Int_t binMax_y = fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->FindBin(+3*fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->GetRMS());
      for (Int_t bin=1; bin<=fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->GetNbinsX(); bin++) {
	if (bin<binMin_y || bin>binMax_y) fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->SetBinContent(bin, 0);
      }
      fHistOffsetSingleMuonsX_vsPtRapidity->SetBinContent(rapBin+1, ptBin+1, fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->GetRMS());
      fHistOffsetSingleMuonsX_vsPtRapidity->SetBinError(rapBin+1, ptBin+1, fHistOffsetSingleMuonsX_tmp[rapBin][ptBin]->GetRMSError());
      fHistOffsetSingleMuonsY_vsPtRapidity->SetBinContent(rapBin+1, ptBin+1, fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->GetRMS());
      fHistOffsetSingleMuonsY_vsPtRapidity->SetBinError(rapBin+1, ptBin+1, fHistOffsetSingleMuonsY_tmp[rapBin][ptBin]->GetRMSError());
    }
  }

  TFile *fileOut = new TFile(Form("%s/%s",fOutputDir.Data(),outputFileName), "recreate");

  printf("Writing output objects to file %s\n", fileOut->GetName());

  fHistOffsetSingleMuonsX -> Write();
  fHistOffsetSingleMuonsY -> Write();
  fHistOffsetSingleMuons  -> Write();
  fHistWOffsetSingleMuons -> Write();
  fHistErrorSingleMuonsX  -> Write();
  fHistErrorSingleMuonsY  -> Write();

  fHistOffsetSingleMuonsX_vsPtRapidity -> Write();
  fHistOffsetSingleMuonsY_vsPtRapidity -> Write();

//   for (Int_t rapBin=0; rapBin<fNRapBinsOffsetSingleMuons; rapBin++) {
//     for (Int_t ptBin=0; ptBin<fNPtBinsOffsetSingleMuons; ptBin++) {
//       fHistOffsetSingleMuonsX_tmp[rapBin][ptBin] -> Write();
//       fHistOffsetSingleMuonsY_tmp[rapBin][ptBin] -> Write();
//     }
//   }

  fHistSingleMuonsPtRapidity -> Write();

  fHistWOffsetMuonPairs        -> Write();
  fHistMassMuonPairs	       -> Write();
  fHistMassMuonPairsWithoutMFT -> Write();
  fHistMassMuonPairsMC         -> Write();
  fHistRapidityPtMuonPairsMC   -> Write();

  fileOut -> Close();

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BookHistos() {

  fHistOffsetSingleMuonsX = new TH1D("fHistOffsetSingleMuonsX", "Offset for single muons along X",  200, -1000, 1000);
  fHistOffsetSingleMuonsY = new TH1D("fHistOffsetSingleMuonsY", "Offset for single muons along Y",  200, -1000, 1000);
  fHistErrorSingleMuonsX  = new TH1D("fHistErrorSingleMuonsX",  "Coordinate Error for single muons along X",  200, 0, 1000);
  fHistErrorSingleMuonsY  = new TH1D("fHistErrorSingleMuonsY",  "Coordinate Error for single muons along Y",  200, 0, 1000);
  fHistOffsetSingleMuons  = new TH1D("fHistOffsetSingleMuons",  "Offset for single muons",          100, 0, 2000);
  fHistWOffsetSingleMuons = new TH1D("fHistWOffsetSingleMuons", "Weighted Offset for single muons", 300, 0, 15);  

  fHistOffsetSingleMuonsX_vsPtRapidity = new TH2D("fHistOffsetSingleMuonsX_vsPtRapidity", "Offset for single muons along X", 
						  10, -4, -2.5, 10, 0.5, 5.5);
  fHistOffsetSingleMuonsY_vsPtRapidity = new TH2D("fHistOffsetSingleMuonsY_vsPtRapidity", "Offset for single muons along Y", 
						  10, -4, -2.5, 10, 0.5, 5.5);

  for (Int_t rapBin=0; rapBin<fNRapBinsOffsetSingleMuons; rapBin++) {
    for (Int_t ptBin=0; ptBin<fNPtBinsOffsetSingleMuons; ptBin++) {
      fHistOffsetSingleMuonsX_tmp[rapBin][ptBin] = new TH1D(Form("fHistOffsetSingleMuonsX_tmp_%02d_%02d",rapBin,ptBin), "",  1000, -1000, 1000);
      fHistOffsetSingleMuonsY_tmp[rapBin][ptBin] = new TH1D(Form("fHistOffsetSingleMuonsY_tmp_%02d_%02d",rapBin,ptBin), "",  1000, -1000, 1000);
    }
  }

  fHistSingleMuonsPtRapidity = new TH2D("fHistSingleMuonsPtRapidity", "Phase Space for single muons", 10, -4, -2.5, 10, 0.5, 5.5);

  fHistOffsetSingleMuonsX -> SetXTitle("Offset(X)  [#mum]");
  fHistOffsetSingleMuonsY -> SetXTitle("Offset(Y)  [#mum]");
  fHistErrorSingleMuonsX  -> SetXTitle("Err. on track position at z_{vtx} (X)  [#mum]");
  fHistErrorSingleMuonsY  -> SetXTitle("Err. on track position at z_{vtx} (Y)  [#mum]");
  fHistOffsetSingleMuons  -> SetXTitle("Offset  [#mum]");
  fHistWOffsetSingleMuons -> SetXTitle("Weighted Offset");

  fHistOffsetSingleMuonsX_vsPtRapidity -> SetXTitle("y^{#mu}");
  fHistOffsetSingleMuonsY_vsPtRapidity -> SetXTitle("y^{#mu}");
  fHistOffsetSingleMuonsX_vsPtRapidity -> SetYTitle("p_{T}^{#mu}  [GeV/c]");
  fHistOffsetSingleMuonsY_vsPtRapidity -> SetYTitle("p_{T}^{#mu}  [GeV/c]");

  fHistSingleMuonsPtRapidity -> SetXTitle("y^{#mu}");
  fHistSingleMuonsPtRapidity -> SetYTitle("p_{T}^{#mu}  [GeV/c]");

  fHistOffsetSingleMuonsX -> Sumw2();
  fHistOffsetSingleMuonsY -> Sumw2();
  fHistErrorSingleMuonsX  -> Sumw2();
  fHistErrorSingleMuonsY  -> Sumw2();
  fHistOffsetSingleMuons  -> Sumw2();
  fHistWOffsetSingleMuons -> Sumw2();

  fHistOffsetSingleMuonsX_vsPtRapidity -> Sumw2();
  fHistOffsetSingleMuonsY_vsPtRapidity -> Sumw2();
  fHistSingleMuonsPtRapidity           -> Sumw2();
    
  //--------------------------------------------

  fHistWOffsetMuonPairs        = new TH1D("fHistWOffsetMuonPairs",        "Weighted Offset for Muon Pairs", 300, 0, 60);
  fHistMassMuonPairs	       = new TH1D("fHistMassMuonPairs",           "Dimuon Mass (MUON+MFT)",  fNMassBins, fMassMin, fMassMax);
  fHistMassMuonPairsWithoutMFT = new TH1D("fHistMassMuonPairsWithoutMFT", "Dimuon Mass (MUON only)", fNMassBins, fMassMin, fMassMax);
  fHistMassMuonPairsMC         = new TH1D("fHistMassMuonPairsMC",         "Dimuon Mass (MC)",        fNMassBins, fMassMin, fMassMax);
  fHistRapidityPtMuonPairsMC   = new TH2D("fHistRapidityPtMuonPairsMC",   "Dimuon Phase Space (MC)", 100, -4.5, -2., 100, 0., 10.); 

  fHistWOffsetMuonPairs        -> SetXTitle("Weighted Offset");
  fHistMassMuonPairs	       -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsWithoutMFT -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsMC         -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistRapidityPtMuonPairsMC   -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairsMC   -> SetYTitle("p_{T}  [GeV/c]");

  fHistWOffsetMuonPairs        -> Sumw2();
  fHistMassMuonPairs	       -> Sumw2();
  fHistMassMuonPairsWithoutMFT -> Sumw2();
  fHistMassMuonPairsMC         -> Sumw2();
  fHistRapidityPtMuonPairsMC   -> Sumw2();

}

//====================================================================================================================================================

