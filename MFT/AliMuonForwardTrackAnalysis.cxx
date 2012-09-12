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
  fNTracksAnalyzedOfEventAfterCut(0),
  fNPairsAnalyzedOfEventAfterCut(0),
  fHistOffsetSingleMuonsX(0x0),
  fHistOffsetSingleMuonsY(0x0),
  fHistOffsetSingleMuons(0x0),
  fHistWOffsetSingleMuons(0x0),
  fHistErrorSingleMuonsX(0x0),
  fHistErrorSingleMuonsY(0x0),
  fHistZOriginSingleMuonsMC(0x0),
  fHistZROriginSingleMuonsMC(0x0), 
  fHistSingleMuonsPtRapidity(0x0), 
  fHistSingleMuonsOffsetChi2(0x0),
  fHistRapidityPtMuonPairs(0x0),
  fHistMassMuonPairsMCVsPt(0x0),
  fHistMassMuonPairsVsPt(0x0),
  fHistMassMuonPairsWithoutMFTVsPt(0x0),
  fHistMassMuonPairsVsPtLSp(0x0),
  fHistMassMuonPairsWithoutMFTVsPtLSp(0x0),
  fHistMassMuonPairsVsPtLSm(0x0),
  fHistMassMuonPairsWithoutMFTVsPtLSm(0x0),
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt(0x0),
  fHistWOffsetMuonPairsAtPCAVsPt(0x0),
  fHistDistancePrimaryVtxPCAVsPt(0x0),
  fHistPCAQualityVsPt(0x0),
  fHistPseudoProperDecayLengthVsPt(0x0),
  fHistDimuonVtxResolutionXVsPt(0x0),
  fHistDimuonVtxResolutionYVsPt(0x0),
  fHistDimuonVtxResolutionZVsPt(0x0),
  fTrueMass(0.),
  fMassMin(0),
  fMassMax(9999),
  fSingleMuonAnalysis(1),
  fMuonPairAnalysis(1),
  fOption(0),
  fTriggerLevel(0),
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
  fCorrelateCutOnOffsetChi2(kFALSE),
  fMaxChi2SingleMuons(1.e9), 
  fMaxOffsetSingleMuons(1.e9),
  fMaxWOffsetMuonPairsAtPrimaryVtx(1.e9), 
  fMaxWOffsetMuonPairsAtPCA(1.e9), 
  fMaxDistancePrimaryVtxPCA(1.e9), 
  fMinPCAQuality(0.)
{

  // default constructor

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::Init(Char_t *inputFileName) {

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

  if (fFirstEvent>=fInputTreeWithBranson->GetEntriesFast()) return kFALSE;
  else if (fFirstEvent<0 || fLastEvent<0 || fFirstEvent>fLastEvent || fFirstEvent>=fInputTreeWithBranson->GetEntriesFast()) {
    fFirstEvent = 0;
    fLastEvent  = fInputTreeWithBranson->GetEntriesFast()-1;
  }
  else {
    fLastEvent = TMath::Min(fLastEvent, Int_t(fInputTreeWithBranson->GetEntriesFast()-1));
  }

  AliInfo(Form("Analysing events %d to %d", fFirstEvent, fLastEvent));

  fMuonForwardTracksWithBranson = new TClonesArray("AliMuonForwardTrack",30);
  fInputTreeWithBranson->SetBranchAddress("tracks", &fMuonForwardTracksWithBranson);  

  fMuonForwardTracksWithoutBranson = new TClonesArray("AliMuonForwardTrack",30);
  fInputTreeWithoutBranson->SetBranchAddress("tracks", &fMuonForwardTracksWithoutBranson);  

  TGeoManager::Import(Form("%s/geometry.root",fInputDir.Data()));

  AliMUONTrackExtrap::SetField();

  fMuonForwardTrackPairsWithBranson    = new TClonesArray("AliMuonForwardTrackPair",10);
  fMuonForwardTrackPairsWithoutBranson = new TClonesArray("AliMuonForwardTrackPair",10);
  fMuonForwardTrackPairsWithBranson    -> SetOwner(kTRUE);
  fMuonForwardTrackPairsWithoutBranson -> SetOwner(kTRUE);
  
  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::LoadNextEvent() {

  if (fEv>fLastEvent) return kFALSE;
  if (fEv<fFirstEvent) { fEv++; return kTRUE; }
  fMuonForwardTracksWithBranson -> Clear("");
  fMuonForwardTracksWithoutBranson -> Clear("");
  fInputTreeWithBranson->GetEvent(fEv);
  fInputTreeWithoutBranson->GetEvent(fEv);
  AliDebug(2,Form("**** analyzing event # %4d (%3d tracks) ****", fEv, fMuonForwardTracksWithBranson->GetEntriesFast()));

  AliInfo(Form("**** analyzing event # %6d of %6d ****", fEv, fLastEvent));

  fPrimaryVtxX = gRandom->Gaus(0., fXVertResMC);
  fPrimaryVtxY = gRandom->Gaus(0., fYVertResMC);
  fPrimaryVtxZ = gRandom->Gaus(0., fZVertResMC);

  if (fSingleMuonAnalysis) {
    fNTracksAnalyzedOfEvent = 0;
    fNTracksAnalyzedOfEventAfterCut = 0;
    fNTracksOfEvent = fMuonForwardTracksWithBranson->GetEntriesFast();
    while (AnalyzeSingleMuon()) continue;
  }
  
  if (fMuonPairAnalysis) {
    if (fMuonForwardTrackPairsWithBranson) {
      fMuonForwardTrackPairsWithBranson->Clear("C");
      fMuonForwardTrackPairsWithoutBranson->Clear("C");
    }
    BuildMuonPairs();
    fNPairsAnalyzedOfEvent = 0;
    fNPairsAnalyzedOfEventAfterCut = 0;
    fNPairsOfEvent = fMuonForwardTrackPairsWithBranson->GetEntriesFast();
    while (AnalyzeMuonPair()) continue;
  }

  AliDebug(2,Form("**** analyzed  event # %4d (%3d tracks and %3d pairs analyzed) ****", fEv, fNTracksAnalyzedOfEventAfterCut, fNPairsAnalyzedOfEventAfterCut));

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
  if (fMFTTrack->GetMatchTrigger() < fTriggerLevel) return kTRUE;
  fMCRefTrack = fMFTTrack->GetMCTrackRef();

  if (fMaxNWrongClustersMC<6) {
    if (!fMCRefTrack) return kTRUE;
    if (fMFTTrack->GetNWrongClustersMC()>fMaxNWrongClustersMC) return kTRUE;
  }

  if (fMCRefTrack) {
    fHistZOriginSingleMuonsMC  -> Fill(-1.*fMCRefTrack->Vz());
    fHistZROriginSingleMuonsMC -> Fill(-1.*fMCRefTrack->Vz(), TMath::Sqrt(fMCRefTrack->Vx()*fMCRefTrack->Vx()+fMCRefTrack->Vy()*fMCRefTrack->Vy()));
  }

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

  fHistSingleMuonsPtRapidity -> Fill(pMu.Rapidity(), pMu.Pt());
  fHistOffsetSingleMuons     -> Fill(1.e4*offset);
  fHistWOffsetSingleMuons    -> Fill(weightedOffset);
  Double_t chi2OverNdf = fMFTTrack->GetGlobalChi2()/Double_t(fMFTTrack->GetNMFTClusters()+fMFTTrack->GetNMUONClusters());
  fHistSingleMuonsOffsetChi2  -> Fill(1.e4*offset, chi2OverNdf);

  fNTracksAnalyzed++;
  fNTracksAnalyzedOfEventAfterCut++;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeMuonPair() {

  if (fNPairsAnalyzedOfEvent>=fNPairsOfEvent) return kFALSE;

  fMFTTrackPairWithBranson    = (AliMuonForwardTrackPair*) fMuonForwardTrackPairsWithBranson->At(fNPairsAnalyzedOfEvent);
  fMFTTrackPairWithoutBranson = (AliMuonForwardTrackPair*) fMuonForwardTrackPairsWithoutBranson->At(fNPairsAnalyzedOfEvent);

  fNPairsAnalyzedOfEvent++;

  if (fUseBransonForKinematics) fMFTTrackPair = fMFTTrackPairWithBranson;
  else fMFTTrackPair = fMFTTrackPairWithoutBranson;

  if ( fMFTTrackPair->GetTrack(0)->GetNWrongClustersMC()>fMaxNWrongClustersMC || 
       fMFTTrackPair->GetTrack(1)->GetNWrongClustersMC()>fMaxNWrongClustersMC ) return kTRUE;

  if (fOption==kResonanceOnly && !fMFTTrackPair->IsResonance())   return kTRUE;
  if (fOption==kResonanceOnly && fMFTTrackPair->GetCharge() != 0) return kTRUE;

  Double_t pca[3] = {0};
  fMFTTrackPair -> GetPointOfClosestApproach(pca);
  Double_t distancePrimaryVtxPCA = TMath::Sqrt(TMath::Power(fPrimaryVtxX-pca[0],2)+
					       TMath::Power(fPrimaryVtxY-pca[1],2)+
					       TMath::Power(fPrimaryVtxZ-pca[2],2));

  fMFTTrackPair -> SetKinem(fPrimaryVtxZ);

  Bool_t passedCut = kFALSE;
  if (fUseBransonForCut) passedCut = PassedCutMuonPair(fMFTTrackPairWithBranson);
  else passedCut = PassedCutMuonPair(fMFTTrackPairWithoutBranson);

  if (!passedCut) return kTRUE;

  // --------------- Filling dimuon histograms --------------------------------

  if (fMFTTrackPair->GetCharge() == 0) {
    if (fOption==kResonanceOnly) fHistMassMuonPairsMCVsPt -> Fill(fMFTTrackPair->GetMassMC(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsVsPt                -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPt      -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
    fHistWOffsetMuonPairsAtPCAVsPt        -> Fill(fMFTTrackPair->GetWeightedOffsetAtPCA(), fMFTTrackPair->GetPt());
    fHistDistancePrimaryVtxPCAVsPt        -> Fill(distancePrimaryVtxPCA*1.e4, fMFTTrackPair->GetPt());
    fHistPCAQualityVsPt                   -> Fill(fMFTTrackPair->GetPCAQuality(), fMFTTrackPair->GetPt());
    if (fOption==kResonanceOnly) fHistPseudoProperDecayLengthVsPt->Fill(GetPseudoProperDecayLength(fMFTTrackPair, fTrueMass), fMFTTrackPair->GetPt());
    if (fEvalDimuonVtxResolution) {
      fHistDimuonVtxResolutionXVsPt->Fill(pca[0]*1.e4, fMFTTrackPair->GetPt());
      fHistDimuonVtxResolutionYVsPt->Fill(pca[1]*1.e4, fMFTTrackPair->GetPt());
      fHistDimuonVtxResolutionZVsPt->Fill(pca[2]*1.e4, fMFTTrackPair->GetPt());
    }
    fHistRapidityPtMuonPairs -> Fill(fMFTTrackPair->GetRapidity(), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == -2) {
    fHistMassMuonPairsVsPtLSm           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSm -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == 2) {
    fHistMassMuonPairsVsPtLSp           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSp -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  }
  
  AliDebug(1, Form("mass = %f   MC = %f", fMFTTrackPair->GetMass(), fMFTTrackPair->GetMassMC()));

  fNPairsAnalyzedOfEventAfterCut++;

  return kTRUE;

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BuildMuonPairs() {

  Int_t nMuonPairs = 0;

  for (Int_t iTrack=0; iTrack<fMuonForwardTracksWithBranson->GetEntriesFast(); iTrack++) {
    for (Int_t jTrack=0; jTrack<iTrack; jTrack++) {
      
      AliMuonForwardTrack *track0_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson->At(iTrack);
      AliMuonForwardTrack *track1_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson->At(jTrack);
      
      AliMuonForwardTrack *track0_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson->At(iTrack);
      AliMuonForwardTrack *track1_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson->At(jTrack);

      if (track0_WithBranson->GetMatchTrigger()<fTriggerLevel || track1_WithBranson->GetMatchTrigger()<fTriggerLevel) continue;

      new ((*fMuonForwardTrackPairsWithBranson)[nMuonPairs]) AliMuonForwardTrackPair(track0_WithBranson, track1_WithBranson);
      new ((*fMuonForwardTrackPairsWithoutBranson)[nMuonPairs]) AliMuonForwardTrackPair(track0_WithoutBranson, track1_WithoutBranson);
      
      AliMuonForwardTrackPair *trackPairWithBranson = (AliMuonForwardTrackPair*) fMuonForwardTrackPairsWithBranson->At(nMuonPairs);
      if (!(fOption==kResonanceOnly && !trackPairWithBranson->IsResonance())) nMuonPairs++;
      
    }
  }

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::PassedCutSingleMuon(AliMuonForwardTrack *track) {

  AliMUONTrackParam *param = track->GetTrackParamAtMFTCluster(0);
  AliMUONTrackExtrap::ExtrapToZCov(param, fPrimaryVtxZ);

  if (track->Pt()<fPtMinSingleMuons) return kFALSE;
  
  Double_t offset = 1.e4*track->GetOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);
  Double_t chi2OverNdf = track->GetGlobalChi2() / Double_t(track->GetNMFTClusters()+track->GetNMUONClusters()); 
  offset /= fMaxOffsetSingleMuons;
  chi2OverNdf /= fMaxChi2SingleMuons;

  if (fCorrelateCutOnOffsetChi2) {
    if (TMath::Sqrt(offset*offset + chi2OverNdf*chi2OverNdf) > 1) return kFALSE;
  }
  else {
    if (offset>1 || chi2OverNdf>1) return kFALSE;
  }

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::PassedCutMuonPair(AliMuonForwardTrackPair *pair) {

  if (!PassedCutSingleMuon(pair->GetTrack(0)) || !PassedCutSingleMuon(pair->GetTrack(1))) return kFALSE;

  if (pair->GetMass()>fMassMax || pair->GetMass()<fMassMin) return kFALSE;

  if (pair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ) > fMaxWOffsetMuonPairsAtPrimaryVtx) return kFALSE;
  if (pair->GetWeightedOffsetAtPCA() > fMaxWOffsetMuonPairsAtPCA) return kFALSE;

  Double_t pca[3] = {0};
  pair -> GetPointOfClosestApproach(pca);
  Double_t distancePrimaryVtxPCA = TMath::Sqrt(TMath::Power(fPrimaryVtxX-pca[0],2)+
					       TMath::Power(fPrimaryVtxY-pca[1],2)+
					       TMath::Power(fPrimaryVtxZ-pca[2],2));

  if (distancePrimaryVtxPCA > fMaxDistancePrimaryVtxPCA) return kFALSE;
  if (pair->GetPCAQuality() < fMinPCAQuality) return kFALSE;

  return kTRUE;

}

//====================================================================================================================================================

Double_t AliMuonForwardTrackAnalysis::GetPseudoProperDecayLength(AliMuonForwardTrackPair *pair, Double_t trueMass) {

  Double_t pca[3] = {0};
  pair -> GetPointOfClosestApproach(pca);

  TVector2 vecVertexPCA(pca[0]-fPrimaryVtxX, pca[1]-fPrimaryVtxY);
  TVector2 dimuonPt(pair->GetPx(), pair->GetPy());
  dimuonPt.Unit();
  
  Double_t l_xy = vecVertexPCA*dimuonPt;
  Double_t pseudoProperDecayLength = (l_xy * trueMass / pair->GetPt()) * 1.e4 ; // in micron

  return pseudoProperDecayLength;
  
}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::Terminate(Char_t *outputFileName) {

  TFile *fileOut = new TFile(Form("%s/%s",fOutputDir.Data(),outputFileName), "recreate");

  printf("Writing output objects to file %s\n", fileOut->GetName());

  fHistOffsetSingleMuonsX  -> Write();
  fHistOffsetSingleMuonsY  -> Write();
  fHistErrorSingleMuonsX   -> Write();
  fHistErrorSingleMuonsY   -> Write();
  fHistOffsetSingleMuons   -> Write();
  fHistWOffsetSingleMuons  -> Write();

  fHistSingleMuonsPtRapidity -> Write();
  fHistSingleMuonsOffsetChi2 -> Write();
  fHistZOriginSingleMuonsMC  -> Write();
  fHistZROriginSingleMuonsMC -> Write();

  fHistMassMuonPairsVsPt                -> Write();
  fHistMassMuonPairsWithoutMFTVsPt      -> Write();
  fHistMassMuonPairsVsPtLSp             -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSp   -> Write();
  fHistMassMuonPairsVsPtLSm             -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSm   -> Write();
  fHistMassMuonPairsMCVsPt              -> Write();
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt -> Write();
  fHistWOffsetMuonPairsAtPCAVsPt        -> Write();
  fHistDistancePrimaryVtxPCAVsPt        -> Write();
  fHistPCAQualityVsPt                   -> Write();
  if (fOption==kResonanceOnly) fHistPseudoProperDecayLengthVsPt->Write();
  if (fEvalDimuonVtxResolution) {
    fHistDimuonVtxResolutionXVsPt -> Write();
    fHistDimuonVtxResolutionYVsPt -> Write();
    fHistDimuonVtxResolutionZVsPt -> Write();
  }
  fHistRapidityPtMuonPairs -> Write();

  fileOut -> Close();

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BookHistos() {

  // -------------------- single muons

  fHistOffsetSingleMuonsX = new TH1D("fHistOffsetSingleMuonsX", "Offset for single muons along X",  200, -1000, 1000);
  fHistOffsetSingleMuonsY = new TH1D("fHistOffsetSingleMuonsY", "Offset for single muons along Y",  200, -1000, 1000);
  fHistErrorSingleMuonsX  = new TH1D("fHistErrorSingleMuonsX",  "Coordinate Error for single muons along X",  200, 0, 1000);
  fHistErrorSingleMuonsY  = new TH1D("fHistErrorSingleMuonsY",  "Coordinate Error for single muons along Y",  200, 0, 1000);
  fHistOffsetSingleMuons  = new TH1D("fHistOffsetSingleMuons",  "Offset for single muons",          200, 0, 2000);
  fHistWOffsetSingleMuons = new TH1D("fHistWOffsetSingleMuons", "Weighted Offset for single muons", 300, 0, 15);  

  fHistSingleMuonsPtRapidity = new TH2D("fHistSingleMuonsPtRapidity", "Phase Space for single muons", 100, -4.5, -2., 100, 0., 10.);
  fHistSingleMuonsOffsetChi2 = new TH2D("fHistSingleMuonsOffsetChi2", "Offset vs #chi^{2}/ndf for single muons", 400, 0, 4000, 200, 0, 10);
  fHistZOriginSingleMuonsMC  = new TH1D("fHistZOriginSingleMuonsMC",  "Z origin for single muons (from MC)",   1000, 0., 500.);
  fHistZROriginSingleMuonsMC = new TH2D("fHistZROriginSingleMuonsMC", "Z-R origin for single muons (from MC)", 1000, 0., 500., 1000, 0., 100.);

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

  fHistZOriginSingleMuonsMC  -> SetXTitle("Z  [cm]");
  fHistZROriginSingleMuonsMC -> SetXTitle("Z  [cm]");
  fHistZROriginSingleMuonsMC -> SetXTitle("R  [cm]");

  fHistOffsetSingleMuonsX -> Sumw2();
  fHistOffsetSingleMuonsY -> Sumw2();
  fHistErrorSingleMuonsX  -> Sumw2();
  fHistErrorSingleMuonsY  -> Sumw2();
  fHistOffsetSingleMuons  -> Sumw2();
  fHistWOffsetSingleMuons -> Sumw2();

  fHistSingleMuonsPtRapidity -> Sumw2();
  fHistSingleMuonsOffsetChi2 -> Sumw2();
  fHistZOriginSingleMuonsMC  -> Sumw2();
  fHistZROriginSingleMuonsMC -> Sumw2();

  // -------------------- muon pairs

  Int_t nBinsPt = 20; Double_t ptMin = 0., ptMax = 10.;   // dimuon pt

  fHistMassMuonPairsVsPt	        = new TH2D("fHistMassMuonPairsVsPt", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPt      = new TH2D("fHistMassMuonPairsWithoutMFTVsPt", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSp	        = new TH2D("fHistMassMuonPairsVsPtLSp", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSp   = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSp", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSm	        = new TH2D("fHistMassMuonPairsVsPtLSm", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSm   = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSm", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsMCVsPt              = new TH2D("fHistMassMuonPairsMCVsPt", "Dimuon Mass (MC) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt = new TH2D("fHistWOffsetMuonPairsAtPrimaryVtxVsPt", "Weighted Offset for Muon Pairs at Primary Vertex vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistWOffsetMuonPairsAtPCAVsPt        = new TH2D("fHistWOffsetMuonPairsAtPCAVsPt", "Weighted Offset for Muon Pairs at PCA vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistDistancePrimaryVtxPCAVsPt        = new TH2D("fHistDistancePrimaryVtxPCA_%d", "Distance between PCA and primary vertex vs p_{T}^{#mu#mu}", 1000, 0, 50000, nBinsPt, ptMin, ptMax);
  fHistPCAQualityVsPt                   = new TH2D("fHistPCAQualityVsPt", "PCA Quality vs p_{T}^{#mu#mu}", 200, 0, 1, nBinsPt, ptMin, ptMax);
  fHistPseudoProperDecayLengthVsPt      = new TH2D("fHistPseudoProperDecayLengthVsPt", "Pseudo proper decay length vs p_{T}^{#mu#mu}", 1000, -5000, 5000, nBinsPt, ptMin, ptMax);
  fHistDimuonVtxResolutionXVsPt = new TH2D("fHistDimuonVtxResolutionXVsPt", "Dimuon Vtx offset along X w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
  fHistDimuonVtxResolutionYVsPt = new TH2D("fHistDimuonVtxResolutionYVsPt", "Dimuon Vtx offset along Y w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
  fHistDimuonVtxResolutionZVsPt = new TH2D("fHistDimuonVtxResolutionZVsPt", "Dimuon Vtx offset along Z w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
    
  fHistMassMuonPairsVsPt                -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsWithoutMFTVsPt      -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsVsPtLSp             -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsWithoutMFTVsPtLSp   -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsVsPtLSm             -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsWithoutMFTVsPtLSm   -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsMCVsPt              -> SetXTitle("Mass  [GeV/c^{2}]");    
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt -> SetXTitle("Weighted Offset at Primary Vtx");
  fHistWOffsetMuonPairsAtPCAVsPt        -> SetXTitle("Weighted Offset at PCA");
  fHistDistancePrimaryVtxPCAVsPt        -> SetXTitle("PCA - Primary Vtx [#mum]");
  fHistPCAQualityVsPt                   -> SetXTitle("PCA Quality");
  fHistPseudoProperDecayLengthVsPt      -> SetXTitle("l_{xy}  [#mum]");
  fHistDimuonVtxResolutionXVsPt         -> SetXTitle("Offset  [#mum]");
  fHistDimuonVtxResolutionYVsPt         -> SetXTitle("Offset  [#mum]");
  fHistDimuonVtxResolutionZVsPt         -> SetXTitle("Offset  [#mum]");

  fHistMassMuonPairsVsPt                -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsWithoutMFTVsPt      -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsVsPtLSp             -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsWithoutMFTVsPtLSp   -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsVsPtLSm             -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsWithoutMFTVsPtLSm   -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsMCVsPt              -> SetYTitle("p_{T}  [GeV/c]");
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt -> SetYTitle("p_{T}  [GeV/c]");
  fHistWOffsetMuonPairsAtPCAVsPt        -> SetYTitle("p_{T}  [GeV/c]");
  fHistDistancePrimaryVtxPCAVsPt        -> SetYTitle("p_{T}  [GeV/c]");  
  fHistPCAQualityVsPt                   -> SetYTitle("p_{T}  [GeV/c]");
  fHistPseudoProperDecayLengthVsPt      -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionXVsPt         -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionYVsPt         -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionZVsPt         -> SetYTitle("p_{T}  [GeV/c]");

  fHistMassMuonPairsVsPt                -> Sumw2();
  fHistMassMuonPairsWithoutMFTVsPt      -> Sumw2();
  fHistMassMuonPairsVsPtLSp             -> Sumw2();
  fHistMassMuonPairsWithoutMFTVsPtLSp   -> Sumw2();
  fHistMassMuonPairsVsPtLSm             -> Sumw2();
  fHistMassMuonPairsWithoutMFTVsPtLSm   -> Sumw2();
  fHistMassMuonPairsMCVsPt              -> Sumw2();    
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt -> Sumw2();
  fHistWOffsetMuonPairsAtPCAVsPt        -> Sumw2();
  fHistDistancePrimaryVtxPCAVsPt        -> Sumw2();
  fHistPCAQualityVsPt                   -> Sumw2();
  fHistPseudoProperDecayLengthVsPt      -> Sumw2();
  fHistDimuonVtxResolutionXVsPt         -> Sumw2();
  fHistDimuonVtxResolutionYVsPt         -> Sumw2();
  fHistDimuonVtxResolutionZVsPt         -> Sumw2();

  fHistRapidityPtMuonPairs = new TH2D("fHistRapidityPtMuonPairs", "Dimuon Phase Space (rec)", 20, -4.5, -2., 20, 0., 10.); 
  fHistRapidityPtMuonPairs -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairs -> SetYTitle("p_{T}  [GeV/c]");
  fHistRapidityPtMuonPairs -> Sumw2();

}

//====================================================================================================================================================

