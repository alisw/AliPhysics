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
  fHistSingleMuonsPtRapidityMC(0x0),
  fHistSingleMuonsOffsetChi2(0x0),
  fHistRapidityPtMuonPairs(0x0),
  fHistMassMuonPairsVsPt(0),
  fHistMassMuonPairsWithoutMFTVsPt(0),
  fHistMassMuonPairsVsPtLSp(0),
  fHistMassMuonPairsWithoutMFTVsPtLSp(0),
  fHistMassMuonPairsVsPtLSm(0),
  fHistMassMuonPairsWithoutMFTVsPtLSm(0),
  fEvalDimuonVtxResolution(kFALSE),
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
  fCenterChi2(0.), 
  fScaleOffset(250.), 
  fScaleChi2(1.5),
  fRadiusCut(1.)
{

  // default constructor

  for (Int_t iPtBin=0; iPtBin<fNMaxPtBinsDimuons+1; iPtBin++) {
    fHistMassMuonPairs[iPtBin]                                = NULL;
    fHistMassMuonPairsWithoutMFT[iPtBin]                      = NULL;
    fHistMassMuonPairsMC[iPtBin]                              = NULL;
    fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin]                 = NULL;
    fHistWOffsetMuonPairsAtPCA[iPtBin]                        = NULL;
    fHistDistancePrimaryVtxPCA[iPtBin]                        = NULL;
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] = NULL;
    fHistDimuonVtxResolutionX[iPtBin]                         = NULL;
    fHistDimuonVtxResolutionY[iPtBin]                         = NULL;
    fHistDimuonVtxResolutionZ[iPtBin]                         = NULL;
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

  if (fFirstEvent<0 || fLastEvent<0 || fFirstEvent>fLastEvent || fFirstEvent>=fInputTreeWithBranson->GetEntriesFast()) {
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
	fMuonForwardTrackPairsWithBranson->SetOwner(kTRUE);
	fMuonForwardTrackPairsWithoutBranson->SetOwner(kTRUE);


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
   // fNPairsOfEvent = fMuonForwardTrackPairsWithBranson->GetEntriesFast();
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
  if (fMatchTrigger && !fMFTTrack->GetMatchTrigger()) return kTRUE;
  fMCRefTrack = fMFTTrack->GetMCTrackRef();

  if (fOption!=kPionsKaons && !fMCRefTrack) return kTRUE;
  if (fOption!=kPionsKaons && fMFTTrack->GetNWrongClustersMC()>fMaxNWrongClustersMC) return kTRUE;

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

  if (fOption!=kPionsKaons) fHistSingleMuonsPtRapidityMC-> Fill(fMCRefTrack->Y(), fMCRefTrack->Pt());
  fHistOffsetSingleMuons      -> Fill(1.e4*offset);
  fHistWOffsetSingleMuons     -> Fill(weightedOffset);
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

  Bool_t passedCut = kFALSE;
  if (fUseBransonForCut) passedCut = PassedCutMuonPair(fMFTTrackPairWithBranson);
  else passedCut = PassedCutMuonPair(fMFTTrackPairWithoutBranson);

  if (!passedCut) return kTRUE;

  if (fUseBransonForKinematics) fMFTTrackPair = fMFTTrackPairWithBranson;
  else fMFTTrackPair = fMFTTrackPairWithoutBranson;

  if ( fMFTTrackPair->GetTrack(0)->GetNWrongClustersMC()>fMaxNWrongClustersMC || 
       fMFTTrackPair->GetTrack(1)->GetNWrongClustersMC()>fMaxNWrongClustersMC ) return kTRUE;

  if (fOption==kResonanceOnly && !fMFTTrackPair->IsResonance()) return kTRUE;

  if (fOption==kResonanceOnly && fMFTTrackPair->GetCharge() != 0) return kTRUE;

  Double_t pca[3] = {0};
  fMFTTrackPair -> GetPointOfClosestApproach(pca);
  Double_t distancePrimaryVtxPCA = TMath::Sqrt(TMath::Power(fPrimaryVtxX-pca[0],2)+
					       TMath::Power(fPrimaryVtxY-pca[1],2)+
					       TMath::Power(fPrimaryVtxZ-pca[2],2));

  fMFTTrackPair -> SetKinem(fPrimaryVtxZ);

  Int_t ptBin = fPtAxisDimuons->FindBin(fMFTTrackPair->GetPt());

  if (1<=ptBin && ptBin<=fNPtDimuBins) {
    fHistMassMuonPairs[ptBin]           -> Fill(fMFTTrackPair->GetMass());
    fHistMassMuonPairsWithoutMFT[ptBin] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
    if (fOption!=kPionsKaons) fHistMassMuonPairsMC[ptBin] -> Fill(fMFTTrackPair->GetMassMC());
    fHistWOffsetMuonPairsAtPrimaryVtx[ptBin] -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
    fHistWOffsetMuonPairsAtPCA[ptBin]        -> Fill(fMFTTrackPair->GetWeightedOffset(pca[0], pca[1], pca[2]));
    fHistDistancePrimaryVtxPCA[ptBin]        -> Fill(distancePrimaryVtxPCA*1.e4);
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[ptBin] -> Fill(fMFTTrackPair->GetWeightedOffset(pca[0], pca[1], pca[2]), distancePrimaryVtxPCA*1.e4);
    if (fEvalDimuonVtxResolution) {
      fHistDimuonVtxResolutionX[ptBin]->Fill(pca[0]*1.e4);
      fHistDimuonVtxResolutionY[ptBin]->Fill(pca[1]*1.e4);
      fHistDimuonVtxResolutionZ[ptBin]->Fill(pca[2]*1.e4);
    }
  }
  fHistMassMuonPairs[0]           -> Fill(fMFTTrackPair->GetMass());
  fHistMassMuonPairsWithoutMFT[0] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
  if (fOption!=kPionsKaons) fHistMassMuonPairsMC[0] -> Fill(fMFTTrackPair->GetMassMC());
  fHistWOffsetMuonPairsAtPrimaryVtx[0] -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ));
  fHistWOffsetMuonPairsAtPCA[0]        -> Fill(fMFTTrackPair->GetWeightedOffset(pca[0], pca[1], pca[2]));
  fHistDistancePrimaryVtxPCA[0]        -> Fill(distancePrimaryVtxPCA*1.e4);
  fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[0] -> Fill(fMFTTrackPair->GetWeightedOffset(pca[0], pca[1], pca[2]), distancePrimaryVtxPCA*1.e4);
  if (fEvalDimuonVtxResolution) {
    fHistDimuonVtxResolutionX[0]->Fill(pca[0]*1.e4);
    fHistDimuonVtxResolutionY[0]->Fill(pca[1]*1.e4);
    fHistDimuonVtxResolutionZ[0]->Fill(pca[2]*1.e4);
  }

  if (fMFTTrackPair->GetCharge() == 0) {
    fHistMassMuonPairsVsPt           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPt -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == -2) {
    fHistMassMuonPairsVsPtLSm           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSm -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == 2) {
    fHistMassMuonPairsVsPtLSp           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSp -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  }
  
  if (fOption!=kPionsKaons) fHistRapidityPtMuonPairs -> Fill(fMFTTrackPair->GetRapidityMC(), fMFTTrackPair->GetPt());
  else fHistRapidityPtMuonPairs -> Fill(fMFTTrackPair->GetRapidity(), fMFTTrackPair->GetPt()); 

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

      if (fMatchTrigger) if (!track0_WithBranson->GetMatchTrigger() || !track1_WithBranson->GetMatchTrigger()) continue;
      if (fOption!=kPionsKaons && (!track0_WithBranson->GetMCTrackRef() || !track1_WithBranson->GetMCTrackRef())) continue;

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

  fHistSingleMuonsPtRapidityMC -> Write();
  fHistSingleMuonsOffsetChi2 -> Write();

  fHistZOriginSingleMuonsMC  -> Write();
  fHistZROriginSingleMuonsMC -> Write();

  for (Int_t iPtBin=0; iPtBin<fNPtDimuBins+1; iPtBin++) {
    fHistMassMuonPairs[iPtBin]                                -> Write();
    fHistMassMuonPairsWithoutMFT[iPtBin]                      -> Write();
    fHistMassMuonPairsMC[iPtBin]                              -> Write();
    fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin]                 -> Write();
    fHistWOffsetMuonPairsAtPCA[iPtBin]                        -> Write();
    fHistDistancePrimaryVtxPCA[iPtBin]                        -> Write();
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] -> Write();
    if (fEvalDimuonVtxResolution) {
      fHistDimuonVtxResolutionX[iPtBin] -> Write();
      fHistDimuonVtxResolutionY[iPtBin] -> Write();
      fHistDimuonVtxResolutionZ[iPtBin] -> Write();
    }
  }

  fHistMassMuonPairsVsPt           -> Write();
  fHistMassMuonPairsWithoutMFTVsPt -> Write();
  fHistMassMuonPairsVsPtLSp           -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSp -> Write();
  fHistMassMuonPairsVsPtLSm           -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSm -> Write();
	
  fHistRapidityPtMuonPairs -> Write();

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

  fHistSingleMuonsPtRapidityMC = new TH2D("fHistSingleMuonsPtRapidityMC", "Phase Space for single muons", 100, -4.5, -2., 100, 0., 10.);
  fHistSingleMuonsOffsetChi2 = new TH2D("fHistSingleMuonsOffsetChi2", "Offset vs #chi^{2}/ndf for single muons", 400, 0, 4000, 200, 0, 10);
  fHistZOriginSingleMuonsMC  = new TH1D("fHistZOriginSingleMuonsMC",  "Z origin for single muons (from MC)",   1000, 0., 500.);
  fHistZROriginSingleMuonsMC = new TH2D("fHistZROriginSingleMuonsMC", "Z-R origin for single muons (from MC)", 1000, 0., 500., 1000, 0., 100.);

  fHistOffsetSingleMuonsX -> SetXTitle("Offset(X)  [#mum]");
  fHistOffsetSingleMuonsY -> SetXTitle("Offset(Y)  [#mum]");
  fHistErrorSingleMuonsX  -> SetXTitle("Err. on track position at z_{vtx} (X)  [#mum]");
  fHistErrorSingleMuonsY  -> SetXTitle("Err. on track position at z_{vtx} (Y)  [#mum]");
  fHistOffsetSingleMuons  -> SetXTitle("Offset  [#mum]");
  fHistWOffsetSingleMuons -> SetXTitle("Weighted Offset");

  fHistSingleMuonsPtRapidityMC -> SetXTitle("y^{#mu}");
  fHistSingleMuonsPtRapidityMC -> SetYTitle("p_{T}^{#mu}  [GeV/c]");
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

  fHistZOriginSingleMuonsMC  -> Sumw2();
  fHistZROriginSingleMuonsMC -> Sumw2();

  fHistSingleMuonsPtRapidityMC -> Sumw2();
  fHistSingleMuonsOffsetChi2 -> Sumw2();
    
  //--------------------------------------------

  for (Int_t iPtBin=0; iPtBin<=fNPtDimuBins+1; iPtBin++) {

    if (!iPtBin) {
      fHistMassMuonPairs[iPtBin]	   = new TH1D(Form("fHistMassMuonPairs_%d",iPtBin),           
						      "Dimuon Mass (MUON+MFT) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsWithoutMFT[iPtBin] = new TH1D(Form("fHistMassMuonPairsWithoutMFT_%d",iPtBin), 
						      "Dimuon Mass (MUON only) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsMC[iPtBin]         = new TH1D(Form("fHistMassMuonPairsMC_%d",iPtBin),         
						      "Dimuon Mass (MC) (All p_{T}^{#mu#mu})",
						      fNMassBins, fMassMin, fMassMax);
      fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin] = new TH1D(Form("fHistWOffsetMuonPairsAtPrimaryVtx_%d",iPtBin),        
							   "Weighted Offset for Muon Pairs at Primary Vertex (All p_{T}^{#mu#mu})",
							   300, 0, 60);
      fHistWOffsetMuonPairsAtPCA[iPtBin]        = new TH1D(Form("fHistWOffsetMuonPairsAtPCA_%d",iPtBin),        
							   "Weighted Offset for Muon Pairs at PCA (All p_{T}^{#mu#mu})",
							   300, 0, 60);
      fHistDistancePrimaryVtxPCA[iPtBin]        = new TH1D(Form("fHistDistancePrimaryVtxPCA_%d",iPtBin),        
							   "Distance between PCA and primary vertex (All p_{T}^{#mu#mu})",
							   1000, 0, 50000);
      fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] = new TH2D(Form("fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA_%d",iPtBin),
									   "DistancePrimaryVtxPCA vs WOffsetMuonPairsAtPCA (All p_{T}^{#mu#mu})",
									   300, 0, 60, 1000, 0, 50000);
      if (fEvalDimuonVtxResolution) {
	fHistDimuonVtxResolutionX[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionX_%d",iPtBin),        
						     "Dimuon Vtx offset along X (All p_{T}^{#mu#mu})",
						     2000, -1000, 1000);
	fHistDimuonVtxResolutionY[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionY_%d",iPtBin),        
						     "Dimuon Vtx offset along Y (All p_{T}^{#mu#mu})",
						     2000, -1000, 1000);
	fHistDimuonVtxResolutionZ[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionZ_%d",iPtBin),        
						     "Dimuon Vtx offset along Z (All p_{T}^{#mu#mu})",
						     2000, -10000, 10000);
      }
    }

    else {
      Double_t ptMin = fPtAxisDimuons->GetBinLowEdge(iPtBin);
      Double_t ptMax = fPtAxisDimuons->GetBinUpEdge(iPtBin);
      fHistMassMuonPairs[iPtBin]	   = new TH1D(Form("fHistMassMuonPairs_%d",iPtBin),           
						      Form("Dimuon Mass (MUON+MFT) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsWithoutMFT[iPtBin] = new TH1D(Form("fHistMassMuonPairsWithoutMFT_%d",iPtBin), 
						      Form("Dimuon Mass (MUON only) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);
      fHistMassMuonPairsMC[iPtBin]         = new TH1D(Form("fHistMassMuonPairsMC_%d",iPtBin),         
						      Form("Dimuon Mass (MC) (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax), 
						      fNMassBins, fMassMin, fMassMax);
      fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin] = new TH1D(Form("fHistWOffsetMuonPairsAtPrimaryVtx_%d",iPtBin),        
							   Form("Weighted Offset for Muon Pairs at Primary Vertex (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",
								ptMin,ptMax), 300, 0, 60);
      fHistWOffsetMuonPairsAtPCA[iPtBin]        = new TH1D(Form("fHistWOffsetMuonPairsAtPCA_%d",iPtBin),        
							   Form("Weighted Offset for Muon Pairs at PCA (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",
								ptMin,ptMax), 300, 0, 60);
      fHistDistancePrimaryVtxPCA[iPtBin]        = new TH1D(Form("fHistDistancePrimaryVtxPCA_%d",iPtBin),        
							   Form("Distance between PCA and primary vertex (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",
								ptMin,ptMax), 1000, 0, 50000);
      fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] = new TH2D(Form("fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA_%d",iPtBin),
									   Form("DistancePrimaryVtxPCA vs WOffsetMuonPairsAtPCA (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",
										ptMin,ptMax), 300, 0, 60, 1000, 0, 50000);
      if (fEvalDimuonVtxResolution) {
	fHistDimuonVtxResolutionX[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionX_%d",iPtBin),        
						     Form("Dimuon Vtx offset along X (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax),
						     2000, -1000, 1000);
	fHistDimuonVtxResolutionY[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionY_%d",iPtBin),        
						     Form("Dimuon Vtx offset along Y (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax),
						     2000, -1000, 1000);
	fHistDimuonVtxResolutionZ[iPtBin] = new TH1D(Form("fHistDimuonVtxResolutionZ_%d",iPtBin),        
						     Form("Dimuon Vtx offset along Z (%3.1f < p_{T}^{#mu#mu} < %3.1f GeV/c)",ptMin,ptMax),
						     2000, -10000, 10000);
      }
    }
    
    fHistMassMuonPairs[iPtBin]                -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsWithoutMFT[iPtBin]      -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsMC[iPtBin]              -> SetXTitle("Mass  [GeV/c^{2}]");    
    fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin] -> SetXTitle("Weighted Offset at Primary Vtx");
    fHistWOffsetMuonPairsAtPCA[iPtBin]        -> SetXTitle("Weighted Offset at PCA");
    fHistDistancePrimaryVtxPCA[iPtBin]        -> SetXTitle("PCA - Primary Vtx [#mum]");
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] -> SetXTitle("Weighted Offset at PCA");
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] -> SetYTitle("PCA - Primary Vtx [#mum]");
    if (fEvalDimuonVtxResolution) {
      fHistDimuonVtxResolutionX[iPtBin] -> SetXTitle("Offset  [#mum]");
      fHistDimuonVtxResolutionY[iPtBin] -> SetXTitle("Offset  [#mum]");
      fHistDimuonVtxResolutionZ[iPtBin] -> SetXTitle("Offset  [#mum]");
    }

    fHistMassMuonPairs[iPtBin]                -> Sumw2();
    fHistMassMuonPairsWithoutMFT[iPtBin]      -> Sumw2();
    fHistMassMuonPairsMC[iPtBin]              -> Sumw2();    
    fHistWOffsetMuonPairsAtPrimaryVtx[iPtBin] -> Sumw2();
    fHistWOffsetMuonPairsAtPCA[iPtBin]        -> Sumw2();
    fHistDistancePrimaryVtxPCA[iPtBin]        -> Sumw2();
    fHistDistancePrimaryVtxPCAvsWOffsetMuonPairsAtPCA[iPtBin] -> Sumw2();
    if (fEvalDimuonVtxResolution) {
      fHistDimuonVtxResolutionX[iPtBin] -> Sumw2();
      fHistDimuonVtxResolutionY[iPtBin] -> Sumw2();
      fHistDimuonVtxResolutionZ[iPtBin] -> Sumw2();
    }

  }
  
  fHistMassMuonPairsVsPt           = new TH2D("fHistMassMuonPairsVsPt",           "Dimuon Mass (MUON+MFT) vs p_{T}",  fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
  fHistMassMuonPairsWithoutMFTVsPt = new TH2D("fHistMassMuonPairsWithoutMFTVsPt", "Dimuon Mass (MUON only) vs p_{T}", fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
	fHistMassMuonPairsVsPtLSp           = new TH2D("fHistMassMuonPairsVsPtLSp",           "Dimuon Mass (MUON+MFT) vs p_{T} Like sign ++",  fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
	fHistMassMuonPairsWithoutMFTVsPtLSp = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSp", "Dimuon Mass (MUON only) vs p_{T} Like sign ++", fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
	fHistMassMuonPairsVsPtLSm           = new TH2D("fHistMassMuonPairsVsPtLSm",           "Dimuon Mass (MUON+MFT) vs p_{T} Like sign --",  fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
	fHistMassMuonPairsWithoutMFTVsPtLSm = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSm", "Dimuon Mass (MUON only) vs p_{T} Like sign --", fNMassBins, fMassMin, fMassMax, 20, 0., 10.);
	
  fHistMassMuonPairsVsPt           -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsWithoutMFTVsPt -> SetXTitle("Mass  [GeV/c^{2}]");
  fHistMassMuonPairsVsPt           -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsWithoutMFTVsPt -> SetYTitle("p_{T}  [GeV/c]");
  fHistMassMuonPairsVsPt           -> Sumw2();
  fHistMassMuonPairsWithoutMFTVsPt -> Sumw2();
	fHistMassMuonPairsVsPtLSp           -> SetXTitle("Mass  [GeV/c^{2}]");
	fHistMassMuonPairsWithoutMFTVsPtLSp -> SetXTitle("Mass  [GeV/c^{2}]");
	fHistMassMuonPairsVsPtLSp           -> SetYTitle("p_{T}  [GeV/c]");
	fHistMassMuonPairsWithoutMFTVsPtLSp -> SetYTitle("p_{T}  [GeV/c]");
	fHistMassMuonPairsVsPtLSp           -> Sumw2();
	fHistMassMuonPairsWithoutMFTVsPtLSp -> Sumw2();
	fHistMassMuonPairsVsPtLSm           -> SetXTitle("Mass  [GeV/c^{2}]");
	fHistMassMuonPairsWithoutMFTVsPtLSm -> SetXTitle("Mass  [GeV/c^{2}]");
	fHistMassMuonPairsVsPtLSm           -> SetYTitle("p_{T}  [GeV/c]");
	fHistMassMuonPairsWithoutMFTVsPtLSm -> SetYTitle("p_{T}  [GeV/c]");
	fHistMassMuonPairsVsPtLSm           -> Sumw2();
	fHistMassMuonPairsWithoutMFTVsPtLSm -> Sumw2();
	
  if (fOption==kPionsKaons) fHistRapidityPtMuonPairs = new TH2D("fHistRapidityPtMuonPairs", "Dimuon Phase Space (rec)", 20, -4.5, -2., 20, 0., 10.); 
  else                      fHistRapidityPtMuonPairs = new TH2D("fHistRapidityPtMuonPairs", "Dimuon Phase Space (MC)", 100, -4.5, -2., 100, 0., 10.); 
  fHistRapidityPtMuonPairs   -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairs   -> SetYTitle("p_{T}  [GeV/c]");
  fHistRapidityPtMuonPairs   -> Sumw2();

}

//====================================================================================================================================================

