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
#include "TGeoManager.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TGraph.h"
#include "AliMuonForwardTrackAnalysis.h"
#include "TSystem.h"

ClassImp(AliMuonForwardTrackAnalysis)

//====================================================================================================================================================

AliMuonForwardTrackAnalysis::AliMuonForwardTrackAnalysis():
  TObject(),
  fInputDir(0),
  fOutputDir(0),
  fInputTreeWithBranson(0x0),
  fInputTreeWithoutBranson(0x0),
  fMuonForwardTracksWithBranson(0),
  fMuonForwardTracksWithBransonMix(0),
  fMuonForwardTrackPairsWithBranson(0),
  fMuonForwardTracksWithoutBranson(0),
  fMuonForwardTracksWithoutBransonMix(0),
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
  fHistXOffsetSingleMuonsVsEtaVsP(0x0),
  fHistYOffsetSingleMuonsVsEtaVsP(0x0),
  fHistOffsetSingleMuonsVsEtaVsP(0x0),
  fHistWOffsetSingleMuonsVsEtaVsP(0x0),
  fHistXOffsetSingleMuonsVsEtaVsPt(0x0),
  fHistYOffsetSingleMuonsVsEtaVsPt(0x0),
  fHistOffsetSingleMuonsVsEtaVsPt(0x0),
  fHistWOffsetSingleMuonsVsEtaVsPt(0x0),
  fHistXErrorSingleMuonsVsEtaVsP(0x0),
  fHistYErrorSingleMuonsVsEtaVsP(0x0),
  fHistXErrorSingleMuonsVsEtaVsPt(0x0),
  fHistYErrorSingleMuonsVsEtaVsPt(0x0),
  fHistZOriginSingleMuonsMC(0x0),
  fHistZROriginSingleMuonsMC(0x0), 
  fHistSingleMuonsPtRapidity(0x0), 
  fHistSingleMuonsOffsetChi2(0x0),
  fHistMassMuonPairsMCVsPt(0x0),
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
  fMinPtSingleMuons(0),
  fMinEtaSingleMuons(-99999.),
  fMaxEtaSingleMuons(+99999.),
  fUseBransonForCut(kFALSE),
  fUseBransonForKinematics(kFALSE),
  fCorrelateCutOnOffsetChi2(kFALSE),
  fMaxChi2SingleMuons(1.e9), 
  fMaxOffsetSingleMuons(1.e9),
  fMaxWOffsetMuonPairsAtPrimaryVtx(1.e9), 
  fMaxWOffsetMuonPairsAtPCA(1.e9), 
  fMaxDistancePrimaryVtxPCA(1.e9), 
  fMinPCAQuality(0.),
  fMixing(kFALSE),
  fNEventsToMix(10)
{

  // default constructor

  for (Int_t i=0; i<2; i++) {
    fHistRapidityPtMuonPairs[i]              = 0;
    fHistMassMuonPairsVsPt[i]                = 0;
    fHistMassMuonPairsWithoutMFTVsPt[i]      = 0;
    fHistMassMuonPairsVsPtLSp[i]             = 0;
    fHistMassMuonPairsWithoutMFTVsPtLSp[i]   = 0;
    fHistMassMuonPairsVsPtLSm[i]             = 0;
    fHistMassMuonPairsWithoutMFTVsPtLSm[i]   = 0;
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[i] = 0;
    fHistWOffsetMuonPairsAtPCAVsPt[i]        = 0;
    fHistDistancePrimaryVtxPCAVsPt[i]        = 0;
    fHistPCAQualityVsPt[i]                   = 0;
    fHistPseudoProperDecayLengthVsPt[i]      = 0;
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
  
  fInputTreeWithBranson    = (TTree*) inputFile->Get("AliMuonForwardTracksWithBranson");
  if (!fInputTreeWithBranson) {
    AliError("Error reading input tree");
    return kFALSE;
  }
  
  fInputTreeWithoutBranson    = (TTree*) inputFile->Get("AliMuonForwardTracksWithoutBranson");
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

  fMuonForwardTracksWithBranson    -> Clear("C");
  fMuonForwardTracksWithoutBranson -> Clear("C");
  fMuonForwardTracksWithBranson    -> Delete();
  fMuonForwardTracksWithoutBranson -> Delete();
  fInputTreeWithBranson    -> GetEvent(fEv);
  fInputTreeWithoutBranson -> GetEvent(fEv);

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
    if (fMuonForwardTrackPairsWithBranson || fMuonForwardTrackPairsWithoutBranson) {
      fMuonForwardTrackPairsWithBranson    -> Clear("C");
      fMuonForwardTrackPairsWithoutBranson -> Clear("C");
      fMuonForwardTrackPairsWithBranson    -> Delete();
      fMuonForwardTrackPairsWithoutBranson -> Delete();
    }
    BuildMuonPairs();
    fNPairsAnalyzedOfEvent = 0;
    fNPairsAnalyzedOfEventAfterCut = 0;
    fNPairsOfEvent = fMuonForwardTrackPairsWithBranson->GetEntriesFast();
    while (AnalyzeMuonPair(kSingleEvent)) continue;
  }

  AliDebug(2,Form("**** analyzed  event # %4d (%3d tracks and %3d pairs analyzed) ****", fEv, fNTracksAnalyzedOfEventAfterCut, fNPairsAnalyzedOfEventAfterCut));
  
  if (fMuonPairAnalysis && fMixing) {
    for (fEvMix=fEv+1; fEvMix<=fEv+fNEventsToMix; fEvMix++) {
      if (fEvMix<=fLastEvent) {

	if (fMuonForwardTracksWithBransonMix || fMuonForwardTracksWithoutBransonMix) {
	  fMuonForwardTracksWithBransonMix    -> Delete();
	  fMuonForwardTracksWithoutBransonMix -> Delete();
	  delete fMuonForwardTracksWithBransonMix;
	  delete fMuonForwardTracksWithoutBransonMix;
	}

	fMuonForwardTracksWithBranson    -> Clear("C");
	fMuonForwardTracksWithoutBranson -> Clear("C");
	fMuonForwardTracksWithBranson    -> Delete();
	fMuonForwardTracksWithoutBranson -> Delete();
	fInputTreeWithBranson    -> GetEvent(fEvMix);
	fInputTreeWithoutBranson -> GetEvent(fEvMix);

	fMuonForwardTracksWithBransonMix    = new TClonesArray(*fMuonForwardTracksWithBranson);
	fMuonForwardTracksWithoutBransonMix = new TClonesArray(*fMuonForwardTracksWithoutBranson);

	fMuonForwardTracksWithBranson    -> Clear("C");
	fMuonForwardTracksWithoutBranson -> Clear("C");
	fMuonForwardTracksWithBranson    -> Delete();
	fMuonForwardTracksWithoutBranson -> Delete();
	fInputTreeWithBranson    -> GetEvent(fEv);
	fInputTreeWithoutBranson -> GetEvent(fEv);

	AliDebug(2,Form("**** mixing event # %4d (%3d tracks) with event # %4d (%3d tracks) ****", 
			fEv,    fMuonForwardTracksWithBranson->GetEntriesFast(),
			fEvMix, fMuonForwardTracksWithBransonMix ->GetEntriesFast()));	  
	if (fMuonForwardTrackPairsWithBranson || fMuonForwardTrackPairsWithoutBranson) {
	  fMuonForwardTrackPairsWithBranson    -> Clear("C");
	  fMuonForwardTrackPairsWithoutBranson -> Clear("C");
	  fMuonForwardTrackPairsWithBranson    -> Delete();
	  fMuonForwardTrackPairsWithoutBranson -> Delete();
	}
	BuildMuonPairsMix();
	fNPairsAnalyzedOfEvent = 0;
	fNPairsAnalyzedOfEventAfterCut = 0;
	fNPairsOfEvent = fMuonForwardTrackPairsWithBranson->GetEntriesFast();
	while (AnalyzeMuonPair(kMixedEvent)) continue;
	AliDebug(2,Form("**** analyzed  mixed event pair (%4d ; %4d) : %3d pairs analyzed ****", fEv, fEvMix, fNPairsAnalyzedOfEventAfterCut));
	
	// delete fMuonForwardTracksWithBransonMix;
	// delete fMuonForwardTracksWithoutBransonMix;

      }
    }
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

  fMFTTrack -> EvalKinem(fPrimaryVtxZ);

  TMatrixD cov(5,5);
  cov = fMFTTrack->GetParamCovMatrix();

  fHistXErrorSingleMuonsVsEtaVsP  -> Fill(1.e4*TMath::Sqrt(cov(0,0)), fMFTTrack->Eta(), fMFTTrack->P());
  fHistYErrorSingleMuonsVsEtaVsP  -> Fill(1.e4*TMath::Sqrt(cov(2,2)), fMFTTrack->Eta(), fMFTTrack->P());
  fHistXErrorSingleMuonsVsEtaVsPt -> Fill(1.e4*TMath::Sqrt(cov(0,0)), fMFTTrack->Eta(), fMFTTrack->Pt());
  fHistYErrorSingleMuonsVsEtaVsPt -> Fill(1.e4*TMath::Sqrt(cov(2,2)), fMFTTrack->Eta(), fMFTTrack->Pt());

  fHistXErrorSingleMuonsVsEtaVsP  -> Fill(1.e4*TMath::Sqrt(cov(0,0)), fMFTTrack->Eta(), fMFTTrack->P());
  fHistYErrorSingleMuonsVsEtaVsP  -> Fill(1.e4*TMath::Sqrt(cov(2,2)), fMFTTrack->Eta(), fMFTTrack->P());
  fHistXErrorSingleMuonsVsEtaVsPt -> Fill(1.e4*TMath::Sqrt(cov(0,0)), fMFTTrack->Eta(), fMFTTrack->Pt());
  fHistYErrorSingleMuonsVsEtaVsPt -> Fill(1.e4*TMath::Sqrt(cov(2,2)), fMFTTrack->Eta(), fMFTTrack->Pt());

  Double_t dX = fMFTTrack->GetOffsetX(fPrimaryVtxX, fPrimaryVtxZ);
  Double_t dY = fMFTTrack->GetOffsetY(fPrimaryVtxY, fPrimaryVtxZ);
  
  Double_t offset = fMFTTrack->GetOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);
  Double_t weightedOffset = fMFTTrack->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ);

  //  AliDebug(2, Form("pdg code = %d\n", fMCRefTrack->GetPdgCode()));

  fHistXOffsetSingleMuonsVsEtaVsP  -> Fill(1.e4*dX,        fMFTTrack->Eta(), fMFTTrack->P());
  fHistYOffsetSingleMuonsVsEtaVsP  -> Fill(1.e4*dY,        fMFTTrack->Eta(), fMFTTrack->P());
  fHistOffsetSingleMuonsVsEtaVsP   -> Fill(1.e4*offset,    fMFTTrack->Eta(), fMFTTrack->P());
  fHistWOffsetSingleMuonsVsEtaVsP  -> Fill(weightedOffset, fMFTTrack->Eta(), fMFTTrack->P());
  fHistXOffsetSingleMuonsVsEtaVsPt -> Fill(1.e4*dX,        fMFTTrack->Eta(), fMFTTrack->Pt());
  fHistYOffsetSingleMuonsVsEtaVsPt -> Fill(1.e4*dY,        fMFTTrack->Eta(), fMFTTrack->Pt());
  fHistOffsetSingleMuonsVsEtaVsPt  -> Fill(1.e4*offset,    fMFTTrack->Eta(), fMFTTrack->Pt());
  fHistWOffsetSingleMuonsVsEtaVsPt -> Fill(weightedOffset, fMFTTrack->Eta(), fMFTTrack->Pt());

  fHistSingleMuonsPtRapidity -> Fill(fMFTTrack->Rapidity(), fMFTTrack->Pt());
  Double_t chi2OverNdf = fMFTTrack->GetGlobalChi2()/Double_t(fMFTTrack->GetNMFTClusters()+fMFTTrack->GetNMUONClusters());
  fHistSingleMuonsOffsetChi2  -> Fill(1.e4*offset, chi2OverNdf);

  fNTracksAnalyzed++;
  fNTracksAnalyzedOfEventAfterCut++;

  return kTRUE;

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::AnalyzeMuonPair(Int_t opt) {

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
    if (fOption==kResonanceOnly && opt==kSingleEvent) fHistMassMuonPairsMCVsPt -> Fill(fMFTTrackPair->GetMassMC(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsVsPt[opt]                -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPt[opt]      -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[opt] -> Fill(fMFTTrackPair->GetWeightedOffset(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
    fHistWOffsetMuonPairsAtPCAVsPt[opt]        -> Fill(fMFTTrackPair->GetWeightedOffsetAtPCA(), fMFTTrackPair->GetPt());
    fHistDistancePrimaryVtxPCAVsPt[opt]        -> Fill(distancePrimaryVtxPCA*1.e4, fMFTTrackPair->GetPt());
    fHistPCAQualityVsPt[opt]                   -> Fill(fMFTTrackPair->GetPCAQuality(), fMFTTrackPair->GetPt());
    fHistPseudoProperDecayLengthVsPt[opt]      -> Fill(GetPseudoProperDecayLength(fMFTTrackPair, fTrueMass), fMFTTrackPair->GetPt());
    if (fEvalDimuonVtxResolution && opt==kSingleEvent) {
      fHistDimuonVtxResolutionXVsPt->Fill(pca[0]*1.e4, fMFTTrackPair->GetPt());
      fHistDimuonVtxResolutionYVsPt->Fill(pca[1]*1.e4, fMFTTrackPair->GetPt());
      fHistDimuonVtxResolutionZVsPt->Fill(pca[2]*1.e4, fMFTTrackPair->GetPt());
    }
    fHistRapidityPtMuonPairs[opt] -> Fill(fMFTTrackPair->GetRapidity(), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == -2) {
    fHistMassMuonPairsVsPtLSm[opt]           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSm[opt] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
  } 
  else if (fMFTTrackPair->GetCharge() == 2) {
    fHistMassMuonPairsVsPtLSp[opt]           -> Fill(fMFTTrackPair->GetMass(), fMFTTrackPair->GetPt());
    fHistMassMuonPairsWithoutMFTVsPtLSp[opt] -> Fill(fMFTTrackPair->GetMassWithoutMFT(fPrimaryVtxX, fPrimaryVtxY, fPrimaryVtxZ), fMFTTrackPair->GetPt());
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

      // Resonance   = both mu from the same resonance
      // Open Charm  = both mu from charmed mesons
      // Open Beauty = both mu from beauty mesons
      // Background1mu  = at least one mu from background
      // Background2mu  = both mu from background

      AliMuonForwardTrackPair *trackPairWithBranson = new AliMuonForwardTrackPair(track0_WithBranson, track1_WithBranson);

      if (fOption==kResonanceOnly && !trackPairWithBranson->IsResonance()) { delete trackPairWithBranson; continue; }
      if (fOption==kNoResonances  &&  trackPairWithBranson->IsResonance()) { delete trackPairWithBranson; continue; }
      
      if (fOption==kCharmOnly  && (!track0_WithBranson->IsFromCharm()  || !track1_WithBranson->IsFromCharm()))  { delete trackPairWithBranson; continue; }

      if (fOption==kBeautyOnly && (!track0_WithBranson->IsFromBeauty() || !track1_WithBranson->IsFromBeauty())) { delete trackPairWithBranson; continue; }

      if (fOption==kBackground1mu && (!track0_WithBranson->IsFromBackground() && !track1_WithBranson->IsFromBackground())) { delete trackPairWithBranson; continue; }
      if (fOption==kBackground2mu && (!track0_WithBranson->IsFromBackground() || !track1_WithBranson->IsFromBackground())) { delete trackPairWithBranson; continue; }

      delete trackPairWithBranson;
      
      new ((*fMuonForwardTrackPairsWithBranson)[nMuonPairs]) AliMuonForwardTrackPair(track0_WithBranson, track1_WithBranson);
      new ((*fMuonForwardTrackPairsWithoutBranson)[nMuonPairs]) AliMuonForwardTrackPair(track0_WithoutBranson, track1_WithoutBranson);

      nMuonPairs++;
      
    }
  }

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BuildMuonPairsMix() {

  Int_t nMuonPairs = 0;

  for (Int_t iTrack=0; iTrack<fMuonForwardTracksWithBranson->GetEntriesFast(); iTrack++) {
    for (Int_t jTrack=0; jTrack<fMuonForwardTracksWithBransonMix->GetEntriesFast(); jTrack++) {
      
      AliMuonForwardTrack *track0_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBranson -> At(iTrack);
      AliMuonForwardTrack *track1_WithBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithBransonMix  -> At(jTrack);
      
      AliMuonForwardTrack *track0_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBranson -> At(iTrack);
      AliMuonForwardTrack *track1_WithoutBranson = (AliMuonForwardTrack*) fMuonForwardTracksWithoutBransonMix  -> At(jTrack);

      if (track0_WithBranson->GetMatchTrigger()<fTriggerLevel || track1_WithBranson->GetMatchTrigger()<fTriggerLevel) continue;

      // Resonance   = both mu from the same resonance
      // Open Charm  = both mu from charmed mesons
      // Open Beauty = both mu from beauty mesons
      // Background1mu  = at least one mu from background
      // Background2mu  = both mu from background

      AliMuonForwardTrackPair *trackPairWithBranson = new AliMuonForwardTrackPair(track0_WithBranson, track1_WithBranson);

      if (fOption==kResonanceOnly && !trackPairWithBranson->IsResonance()) { delete trackPairWithBranson; continue; }
      if (fOption==kNoResonances  &&  trackPairWithBranson->IsResonance()) { delete trackPairWithBranson; continue; }
      
      if (fOption==kCharmOnly  && (!track0_WithBranson->IsFromCharm()  || !track1_WithBranson->IsFromCharm()))  { delete trackPairWithBranson; continue; }

      if (fOption==kBeautyOnly && (!track0_WithBranson->IsFromBeauty() || !track1_WithBranson->IsFromBeauty())) { delete trackPairWithBranson; continue; }

      if (fOption==kBackground1mu && (!track0_WithBranson->IsFromBackground() && !track1_WithBranson->IsFromBackground())) { delete trackPairWithBranson; continue; }
      if (fOption==kBackground2mu && (!track0_WithBranson->IsFromBackground() || !track1_WithBranson->IsFromBackground())) { delete trackPairWithBranson; continue; }

      delete trackPairWithBranson;

      new ((*fMuonForwardTrackPairsWithBranson)[nMuonPairs])    AliMuonForwardTrackPair(track0_WithBranson,    track1_WithBranson);
      new ((*fMuonForwardTrackPairsWithoutBranson)[nMuonPairs]) AliMuonForwardTrackPair(track0_WithoutBranson, track1_WithoutBranson);
      
      nMuonPairs++;
      
    }
  }

}

//====================================================================================================================================================

Bool_t AliMuonForwardTrackAnalysis::PassedCutSingleMuon(AliMuonForwardTrack *track) {

  track -> EvalKinem(fPrimaryVtxZ);
  if (track->Pt()<fMinPtSingleMuons) return kFALSE;
  if (track->Eta()<fMinEtaSingleMuons || track->Eta()>fMaxEtaSingleMuons) return kFALSE;
  
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
  TVector2 dimuonPtUnit = dimuonPt.Unit();
  
  Double_t l_xy = vecVertexPCA*dimuonPtUnit;
  Double_t pseudoProperDecayLength = (l_xy * trueMass / pair->GetPt()) * 1.e4 ; // in micron

  return pseudoProperDecayLength;
  
}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::Terminate(Char_t *outputFileName) {

  TFile *fileOut = new TFile(Form("%s/%s",fOutputDir.Data(),outputFileName), "recreate");

  printf("Writing output objects to file %s\n", fileOut->GetName());

  // single muons

  fHistXOffsetSingleMuonsVsEtaVsP  -> Write();
  fHistYOffsetSingleMuonsVsEtaVsP  -> Write();
  fHistXOffsetSingleMuonsVsEtaVsPt -> Write();
  fHistYOffsetSingleMuonsVsEtaVsPt -> Write();

  fHistXErrorSingleMuonsVsEtaVsP   -> Write();
  fHistYErrorSingleMuonsVsEtaVsP   -> Write();
  fHistXErrorSingleMuonsVsEtaVsPt  -> Write();
  fHistYErrorSingleMuonsVsEtaVsPt  -> Write();

  fHistOffsetSingleMuonsVsEtaVsP   -> Write();
  fHistOffsetSingleMuonsVsEtaVsPt  -> Write();
  fHistWOffsetSingleMuonsVsEtaVsP  -> Write();
  fHistWOffsetSingleMuonsVsEtaVsPt -> Write();

  fHistSingleMuonsPtRapidity  -> Write();
  fHistSingleMuonsOffsetChi2  -> Write();
  fHistZOriginSingleMuonsMC   -> Write();
  fHistZROriginSingleMuonsMC  -> Write();

  // dimuons

  fHistMassMuonPairsVsPt[kSingleEvent]                -> Write();
  fHistMassMuonPairsWithoutMFTVsPt[kSingleEvent]      -> Write();
  fHistMassMuonPairsVsPtLSp[kSingleEvent]             -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSp[kSingleEvent]   -> Write();
  fHistMassMuonPairsVsPtLSm[kSingleEvent]             -> Write();
  fHistMassMuonPairsWithoutMFTVsPtLSm[kSingleEvent]   -> Write();
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt[kSingleEvent] -> Write();
  fHistWOffsetMuonPairsAtPCAVsPt[kSingleEvent]        -> Write();
  fHistDistancePrimaryVtxPCAVsPt[kSingleEvent]        -> Write();
  fHistPCAQualityVsPt[kSingleEvent]                   -> Write();
  fHistPseudoProperDecayLengthVsPt[kSingleEvent]      -> Write();
  fHistRapidityPtMuonPairs[kSingleEvent]              -> Write();
  fHistMassMuonPairsMCVsPt                            -> Write();
  if (fEvalDimuonVtxResolution) {
    fHistDimuonVtxResolutionXVsPt -> Write();
    fHistDimuonVtxResolutionYVsPt -> Write();
    fHistDimuonVtxResolutionZVsPt -> Write();
  }

  if (fMixing) {
    fHistMassMuonPairsVsPt[kMixedEvent]                -> Write();
    fHistMassMuonPairsWithoutMFTVsPt[kMixedEvent]      -> Write();
    fHistMassMuonPairsVsPtLSp[kMixedEvent]             -> Write();
    fHistMassMuonPairsWithoutMFTVsPtLSp[kMixedEvent]   -> Write();
    fHistMassMuonPairsVsPtLSm[kMixedEvent]             -> Write();
    fHistMassMuonPairsWithoutMFTVsPtLSm[kMixedEvent]   -> Write();
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[kMixedEvent] -> Write();
    fHistWOffsetMuonPairsAtPCAVsPt[kMixedEvent]        -> Write();
    fHistDistancePrimaryVtxPCAVsPt[kMixedEvent]        -> Write();
    fHistPCAQualityVsPt[kMixedEvent]                   -> Write();
    fHistPseudoProperDecayLengthVsPt[kMixedEvent]      -> Write();
    fHistRapidityPtMuonPairs[kMixedEvent]              -> Write();
  }

  fileOut -> Close();

}

//====================================================================================================================================================

void AliMuonForwardTrackAnalysis::BookHistos() {

  // -------------------- single muons

  fHistXOffsetSingleMuonsVsEtaVsP  = new TH3D("fHistXOffsetSingleMuonsVsEtaVsP",  "Offset for single muons along X vs #eta vs P",      200, -1000, 1000, 100, -4.5, -2., 100, 0, 100);
  fHistYOffsetSingleMuonsVsEtaVsP  = new TH3D("fHistYOffsetSingleMuonsVsEtaVsP",  "Offset for single muons along Y vs #eta vs P",      200, -1000, 1000, 100, -4.5, -2., 100, 0, 100);
  fHistXOffsetSingleMuonsVsEtaVsPt = new TH3D("fHistXOffsetSingleMuonsVsEtaVsPt", "Offset for single muons along X vs #eta vs p_{T}",  200, -1000, 1000, 100, -4.5, -2., 100, 0,  10);
  fHistYOffsetSingleMuonsVsEtaVsPt = new TH3D("fHistYOffsetSingleMuonsVsEtaVsPt", "Offset for single muons along Y vs #eta vs p_{T}",  200, -1000, 1000, 100, -4.5, -2., 100, 0,  10);

  fHistXErrorSingleMuonsVsEtaVsP   = new TH3D("fHistXErrorSingleMuonsVsEtaVsP",   "Coordinate Error for single muons along X vs #eta vs P",      200, 0, 1000, 100, -4.5, -2., 100, 0, 100);
  fHistYErrorSingleMuonsVsEtaVsP   = new TH3D("fHistYErrorSingleMuonsVsEtaVsP",   "Coordinate Error for single muons along Y vs #eta vs P",      200, 0, 1000, 100, -4.5, -2., 100, 0, 100);
  fHistXErrorSingleMuonsVsEtaVsPt  = new TH3D("fHistXErrorSingleMuonsVsEtaVsPt",  "Coordinate Error for single muons along X vs #eta vs p_{T}",  200, 0, 1000, 100, -4.5, -2., 100, 0,  10);
  fHistYErrorSingleMuonsVsEtaVsPt  = new TH3D("fHistYErrorSingleMuonsVsEtaVsPt",  "Coordinate Error for single muons along Y vs #eta vs p_{T}",  200, 0, 1000, 100, -4.5, -2., 100, 0,  10);

  fHistOffsetSingleMuonsVsEtaVsP   = new TH3D("fHistOffsetSingleMuonsVsEtaVsP",   "Offset for single muons vs #eta vs P",              200, 0, 2000, 100, -4.5, -2., 100, 0, 100);
  fHistOffsetSingleMuonsVsEtaVsPt  = new TH3D("fHistOffsetSingleMuonsVsEtaVsPt",  "Offset for single muons vs #eta vs p_{T}",          200, 0, 2000, 100, -4.5, -2., 100, 0,  10);
  fHistWOffsetSingleMuonsVsEtaVsP  = new TH3D("fHistWOffsetSingleMuonsVsEtaVsP",  "Weighted Offset for single muons vs #eta vs P",     300, 0,   15, 100, -4.5, -2., 100, 0, 100);  
  fHistWOffsetSingleMuonsVsEtaVsPt = new TH3D("fHistWOffsetSingleMuonsVsEtaVsPt", "Weighted Offset for single muons vs #eta vs p_{T}", 300, 0,   15, 100, -4.5, -2., 100, 0,  10);  

  fHistSingleMuonsPtRapidity = new TH2D("fHistSingleMuonsPtRapidity", "Phase Space for single muons", 100, -4.5, -2., 100, 0., 10.);
  fHistSingleMuonsOffsetChi2 = new TH2D("fHistSingleMuonsOffsetChi2", "Offset vs #chi^{2}/ndf for single muons", 400, 0, 4000, 200, 0, 10);
  fHistZOriginSingleMuonsMC  = new TH1D("fHistZOriginSingleMuonsMC",  "Z origin for single muons (from MC)",   1000, 0., 500.);
  fHistZROriginSingleMuonsMC = new TH2D("fHistZROriginSingleMuonsMC", "Z-R origin for single muons (from MC)", 1000, 0., 500., 1000, 0., 100.);

  //

  fHistXOffsetSingleMuonsVsEtaVsP  -> SetXTitle("Offset(X)  [#mum]");
  fHistYOffsetSingleMuonsVsEtaVsP  -> SetXTitle("Offset(Y)  [#mum]");
  fHistXOffsetSingleMuonsVsEtaVsPt -> SetXTitle("Offset(X)  [#mum]");
  fHistYOffsetSingleMuonsVsEtaVsPt -> SetXTitle("Offset(Y)  [#mum]");

  fHistXErrorSingleMuonsVsEtaVsP   -> SetXTitle("Err. on track position at z_{vtx} (X)  [#mum]");
  fHistYErrorSingleMuonsVsEtaVsP   -> SetXTitle("Err. on track position at z_{vtx} (Y)  [#mum]");
  fHistXErrorSingleMuonsVsEtaVsPt  -> SetXTitle("Err. on track position at z_{vtx} (X)  [#mum]");
  fHistYErrorSingleMuonsVsEtaVsPt  -> SetXTitle("Err. on track position at z_{vtx} (Y)  [#mum]");

  fHistOffsetSingleMuonsVsEtaVsP   -> SetXTitle("Offset  [#mum]");
  fHistOffsetSingleMuonsVsEtaVsPt  -> SetXTitle("Offset  [#mum]");
  fHistWOffsetSingleMuonsVsEtaVsP  -> SetXTitle("Weighted Offset");
  fHistWOffsetSingleMuonsVsEtaVsPt -> SetXTitle("Weighted Offset");    

  fHistXOffsetSingleMuonsVsEtaVsP  -> SetYTitle("#eta");
  fHistYOffsetSingleMuonsVsEtaVsP  -> SetYTitle("#eta");
  fHistXOffsetSingleMuonsVsEtaVsPt -> SetYTitle("#eta");
  fHistYOffsetSingleMuonsVsEtaVsPt -> SetYTitle("#eta");

  fHistXErrorSingleMuonsVsEtaVsP   -> SetYTitle("#eta");
  fHistYErrorSingleMuonsVsEtaVsP   -> SetYTitle("#eta");
  fHistXErrorSingleMuonsVsEtaVsPt  -> SetYTitle("#eta");
  fHistYErrorSingleMuonsVsEtaVsPt  -> SetYTitle("#eta");

  fHistOffsetSingleMuonsVsEtaVsP   -> SetYTitle("#eta");
  fHistOffsetSingleMuonsVsEtaVsPt  -> SetYTitle("#eta");
  fHistWOffsetSingleMuonsVsEtaVsP  -> SetYTitle("#eta");
  fHistWOffsetSingleMuonsVsEtaVsPt -> SetYTitle("#eta");    

  fHistXOffsetSingleMuonsVsEtaVsP  -> SetZTitle("P  [GeV/c]");
  fHistYOffsetSingleMuonsVsEtaVsP  -> SetZTitle("P  [GeV/c]");
  fHistXOffsetSingleMuonsVsEtaVsPt -> SetZTitle("p_{T}  [GeV/c]");
  fHistYOffsetSingleMuonsVsEtaVsPt -> SetZTitle("p_{T}  [GeV/c]");

  fHistXErrorSingleMuonsVsEtaVsP   -> SetZTitle("P  [GeV/c]");
  fHistYErrorSingleMuonsVsEtaVsP   -> SetZTitle("P  [GeV/c]");
  fHistXErrorSingleMuonsVsEtaVsPt  -> SetZTitle("p_{T}  [GeV/c]");
  fHistYErrorSingleMuonsVsEtaVsPt  -> SetZTitle("p_{T}  [GeV/c]");

  fHistOffsetSingleMuonsVsEtaVsP   -> SetZTitle("P  [GeV/c]");
  fHistOffsetSingleMuonsVsEtaVsPt  -> SetZTitle("p_{T}  [GeV/c]");
  fHistWOffsetSingleMuonsVsEtaVsP  -> SetZTitle("P  [GeV/c]");
  fHistWOffsetSingleMuonsVsEtaVsPt -> SetZTitle("p_{T}  [GeV/c]");    

  fHistSingleMuonsPtRapidity -> SetXTitle("y^{#mu}");
  fHistSingleMuonsPtRapidity -> SetYTitle("p_{T}^{#mu}  [GeV/c]");
  fHistSingleMuonsOffsetChi2 -> SetXTitle("Offset  [#mum]");
  fHistSingleMuonsOffsetChi2 -> SetYTitle("#chi^{2}/ndf");

  fHistZOriginSingleMuonsMC  -> SetXTitle("Z  [cm]");
  fHistZROriginSingleMuonsMC -> SetXTitle("Z  [cm]");
  fHistZROriginSingleMuonsMC -> SetYTitle("R  [cm]");

  //

  fHistXOffsetSingleMuonsVsEtaVsP  -> Sumw2();
  fHistYOffsetSingleMuonsVsEtaVsP  -> Sumw2();
  fHistXOffsetSingleMuonsVsEtaVsPt -> Sumw2();
  fHistYOffsetSingleMuonsVsEtaVsPt -> Sumw2();

  fHistXErrorSingleMuonsVsEtaVsP   -> Sumw2();
  fHistYErrorSingleMuonsVsEtaVsP   -> Sumw2();
  fHistXErrorSingleMuonsVsEtaVsPt  -> Sumw2();
  fHistYErrorSingleMuonsVsEtaVsPt  -> Sumw2();

  fHistOffsetSingleMuonsVsEtaVsP   -> Sumw2();
  fHistOffsetSingleMuonsVsEtaVsPt  -> Sumw2();
  fHistWOffsetSingleMuonsVsEtaVsP  -> Sumw2();
  fHistWOffsetSingleMuonsVsEtaVsPt -> Sumw2();

  fHistSingleMuonsPtRapidity  -> Sumw2();
  fHistSingleMuonsOffsetChi2  -> Sumw2();
  fHistZOriginSingleMuonsMC   -> Sumw2();
  fHistZROriginSingleMuonsMC  -> Sumw2();

  // -------------------- muon pairs

  Int_t nBinsPt = 20; Double_t ptMin = 0., ptMax = 10.;   // dimuon pt

  fHistMassMuonPairsVsPt[kSingleEvent]	              = new TH2D("fHistMassMuonPairsVsPt", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPt[kSingleEvent]      = new TH2D("fHistMassMuonPairsWithoutMFTVsPt", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSp[kSingleEvent]	      = new TH2D("fHistMassMuonPairsVsPtLSp", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSp[kSingleEvent]   = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSp", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSm[kSingleEvent]	      = new TH2D("fHistMassMuonPairsVsPtLSm", "Dimuon Mass (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSm[kSingleEvent]   = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSm", "Dimuon Mass (MUON only) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt[kSingleEvent] = new TH2D("fHistWOffsetMuonPairsAtPrimaryVtxVsPt", "Weighted Offset for Muon Pairs at Primary Vertex vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistWOffsetMuonPairsAtPCAVsPt[kSingleEvent]        = new TH2D("fHistWOffsetMuonPairsAtPCAVsPt", "Weighted Offset for Muon Pairs at PCA vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistDistancePrimaryVtxPCAVsPt[kSingleEvent]        = new TH2D("fHistDistancePrimaryVtxPCAVsPt", "Distance between PCA and primary vertex vs p_{T}^{#mu#mu}", 1000, 0, 50000, nBinsPt, ptMin, ptMax);
  fHistPCAQualityVsPt[kSingleEvent]                   = new TH2D("fHistPCAQualityVsPt", "PCA Quality vs p_{T}^{#mu#mu}", 200, 0, 1, nBinsPt, ptMin, ptMax);
  fHistPseudoProperDecayLengthVsPt[kSingleEvent]      = new TH2D("fHistPseudoProperDecayLengthVsPt", "Pseudo proper decay length vs p_{T}^{#mu#mu}", 1000, -5000, 5000, nBinsPt, ptMin, ptMax);

  fHistMassMuonPairsVsPt[kMixedEvent]	              = new TH2D("fHistMassMuonPairsVsPtMix", "Dimuon Mass Mix (MUON+MFT) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPt[kMixedEvent]       = new TH2D("fHistMassMuonPairsWithoutMFTVsPtMix", "Dimuon Mass Mix (MUON only) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSp[kMixedEvent]	      = new TH2D("fHistMassMuonPairsVsPtLSpMix", "Dimuon Mass Mix (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSp[kMixedEvent]    = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSpMix", "Dimuon Mass Mix (MUON only) vs p_{T}^{#mu#mu} Like-Sign ++", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistMassMuonPairsVsPtLSm[kMixedEvent]	      = new TH2D("fHistMassMuonPairsVsPtLSmMix", "Dimuon Mass Mix (MUON+MFT) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax);
  fHistMassMuonPairsWithoutMFTVsPtLSm[kMixedEvent]    = new TH2D("fHistMassMuonPairsWithoutMFTVsPtLSmMix", "Dimuon Mass Mix (MUON only) vs p_{T}^{#mu#mu} Like-Sign --", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistWOffsetMuonPairsAtPrimaryVtxVsPt[kMixedEvent]  = new TH2D("fHistWOffsetMuonPairsAtPrimaryVtxVsPtMix", "Weighted Offset for Muon Pairs Mix at Primary Vertex vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistWOffsetMuonPairsAtPCAVsPt[kMixedEvent]         = new TH2D("fHistWOffsetMuonPairsAtPCAVsPtMix", "Weighted Offset for Muon Pairs Mix at PCA vs p_{T}^{#mu#mu}", 300, 0, 60, nBinsPt, ptMin, ptMax);
  fHistDistancePrimaryVtxPCAVsPt[kMixedEvent]         = new TH2D("fHistDistancePrimaryVtxPCAVsPtMix", "Distance between PCA and primary vertex Mix vs p_{T}^{#mu#mu}", 1000, 0, 50000, nBinsPt, ptMin, ptMax);
  fHistPCAQualityVsPt[kMixedEvent]                    = new TH2D("fHistPCAQualityVsPtMix", "PCA Quality Mix vs p_{T}^{#mu#mu}", 200, 0, 1, nBinsPt, ptMin, ptMax);
  fHistPseudoProperDecayLengthVsPt[kMixedEvent]       = new TH2D("fHistPseudoProperDecayLengthVsPtMix", "Pseudo proper decay length Mix vs p_{T}^{#mu#mu}", 1000, -5000, 5000, nBinsPt, ptMin, ptMax);

  fHistMassMuonPairsMCVsPt      = new TH2D("fHistMassMuonPairsMCVsPt", "Dimuon Mass (MC) vs p_{T}^{#mu#mu}", 1000, 0, 10, nBinsPt, ptMin, ptMax); 
  fHistDimuonVtxResolutionXVsPt = new TH2D("fHistDimuonVtxResolutionXVsPt", "Dimuon Vtx offset along X w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
  fHistDimuonVtxResolutionYVsPt = new TH2D("fHistDimuonVtxResolutionYVsPt", "Dimuon Vtx offset along Y w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
  fHistDimuonVtxResolutionZVsPt = new TH2D("fHistDimuonVtxResolutionZVsPt", "Dimuon Vtx offset along Z w.r.t. true Vtx vs p_{T}^{#mu#mu}", 2000, -1000, 1000, nBinsPt, ptMin, ptMax);
    
  for (Int_t i=0; i<2; i++) {
    fHistMassMuonPairsVsPt[i]                -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsWithoutMFTVsPt[i]      -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsVsPtLSp[i]             -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsWithoutMFTVsPtLSp[i]   -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsVsPtLSm[i]             -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistMassMuonPairsWithoutMFTVsPtLSm[i]   -> SetXTitle("Mass  [GeV/c^{2}]");
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[i] -> SetXTitle("Weighted Offset at Primary Vtx");
    fHistWOffsetMuonPairsAtPCAVsPt[i]        -> SetXTitle("Weighted Offset at PCA");
    fHistDistancePrimaryVtxPCAVsPt[i]        -> SetXTitle("PCA - Primary Vtx [#mum]");
    fHistPCAQualityVsPt[i]                   -> SetXTitle("PCA Quality");
    fHistPseudoProperDecayLengthVsPt[i]      -> SetXTitle("l_{xy}  [#mum]");
  }
  fHistMassMuonPairsMCVsPt              -> SetXTitle("Mass  [GeV/c^{2}]");    
  fHistDimuonVtxResolutionXVsPt         -> SetXTitle("Offset  [#mum]");
  fHistDimuonVtxResolutionYVsPt         -> SetXTitle("Offset  [#mum]");
  fHistDimuonVtxResolutionZVsPt         -> SetXTitle("Offset  [#mum]");

  for (Int_t i=0; i<2; i++) {
    fHistMassMuonPairsVsPt[i]                -> SetYTitle("p_{T}  [GeV/c]");
    fHistMassMuonPairsWithoutMFTVsPt[i]      -> SetYTitle("p_{T}  [GeV/c]");
    fHistMassMuonPairsVsPtLSp[i]             -> SetYTitle("p_{T}  [GeV/c]");
    fHistMassMuonPairsWithoutMFTVsPtLSp[i]   -> SetYTitle("p_{T}  [GeV/c]");
    fHistMassMuonPairsVsPtLSm[i]             -> SetYTitle("p_{T}  [GeV/c]");
    fHistMassMuonPairsWithoutMFTVsPtLSm[i]   -> SetYTitle("p_{T}  [GeV/c]");
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[i] -> SetYTitle("p_{T}  [GeV/c]");
    fHistWOffsetMuonPairsAtPCAVsPt[i]        -> SetYTitle("p_{T}  [GeV/c]");
    fHistDistancePrimaryVtxPCAVsPt[i]        -> SetYTitle("p_{T}  [GeV/c]");  
    fHistPCAQualityVsPt[i]                   -> SetYTitle("p_{T}  [GeV/c]");
    fHistPseudoProperDecayLengthVsPt[i]      -> SetYTitle("p_{T}  [GeV/c]");
  }
  fHistMassMuonPairsMCVsPt              -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionXVsPt         -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionYVsPt         -> SetYTitle("p_{T}  [GeV/c]");
  fHistDimuonVtxResolutionZVsPt         -> SetYTitle("p_{T}  [GeV/c]");

  for (Int_t i=0; i<2; i++) {
    fHistMassMuonPairsVsPt[i]                -> Sumw2();
    fHistMassMuonPairsWithoutMFTVsPt[i]      -> Sumw2();
    fHistMassMuonPairsVsPtLSp[i]             -> Sumw2();
    fHistMassMuonPairsWithoutMFTVsPtLSp[i]   -> Sumw2();
    fHistMassMuonPairsVsPtLSm[i]             -> Sumw2();
    fHistMassMuonPairsWithoutMFTVsPtLSm[i]   -> Sumw2();
    fHistWOffsetMuonPairsAtPrimaryVtxVsPt[i] -> Sumw2();
    fHistWOffsetMuonPairsAtPCAVsPt[i]        -> Sumw2();
    fHistDistancePrimaryVtxPCAVsPt[i]        -> Sumw2();
    fHistPCAQualityVsPt[i]                   -> Sumw2();
    fHistPseudoProperDecayLengthVsPt[i]      -> Sumw2();
  }
  fHistMassMuonPairsMCVsPt              -> Sumw2();    
  fHistDimuonVtxResolutionXVsPt         -> Sumw2();
  fHistDimuonVtxResolutionYVsPt         -> Sumw2();
  fHistDimuonVtxResolutionZVsPt         -> Sumw2();

  fHistRapidityPtMuonPairs[kSingleEvent] = new TH2D("fHistRapidityPtMuonPairs", "Dimuon Phase Space (rec)", 20, -4.5, -2., 20, 0., 10.); 
  fHistRapidityPtMuonPairs[kSingleEvent] -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairs[kSingleEvent] -> SetYTitle("p_{T}  [GeV/c]");
  fHistRapidityPtMuonPairs[kSingleEvent] -> Sumw2();

  fHistRapidityPtMuonPairs[kMixedEvent] = new TH2D("fHistRapidityPtMuonPairsMix", "Dimuon Phase Space (rec)", 20, -4.5, -2., 20, 0., 10.); 
  fHistRapidityPtMuonPairs[kMixedEvent] -> SetXTitle("Rapidity");
  fHistRapidityPtMuonPairs[kMixedEvent] -> SetYTitle("p_{T}  [GeV/c]");
  fHistRapidityPtMuonPairs[kMixedEvent] -> Sumw2();

}

//====================================================================================================================================================

