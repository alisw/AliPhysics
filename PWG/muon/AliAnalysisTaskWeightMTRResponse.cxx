/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAnalysisTaskWeightMTRResponse.h"

#include "Riostream.h"
#include <vector>

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "THashList.h"
#include "TLorentzVector.h"

// STEER includes
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliInputEventHandler.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// PWG includes
#include "AliAnalysisMuonUtility.h"
#include "AliMergeableCollection.h"
#include "AliUtilityMuonAncestor.h"
#include "AliUtilityDimuonSource.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskWeightMTRResponse) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskWeightMTRResponse::AliAnalysisTaskWeightMTRResponse() :
  AliAnalysisTaskSE(),
  fUtilityMuonAncestor(),
  fUtilityDimuonSource(),
  fMergeableCollection(0x0),
  fResponseType(0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskWeightMTRResponse::AliAnalysisTaskWeightMTRResponse ( const char *name, AliMTRParameterizedResponse* response, Int_t responseType ) :
  AliAnalysisTaskSE(name),
  fResponse(*response),
  fUtilityMuonAncestor(),
  fUtilityDimuonSource(),
  fMergeableCollection(0x0),
  fResponseType(responseType)
{
  //
  /// Constructor
  //
  std::cout << "Response type: ";
  switch ( fResponseType ) {
    case AliMTRParameterizedResponse::kAptOverAll:
    std::cout << "Apt / all reco tracks";
    break;
    case AliMTRParameterizedResponse::kLptOverApt:
    std::cout << "Lpt / Apt";
    break;
    case AliMTRParameterizedResponse::kHptOverLpt:
    std::cout << "Hpt / Lpt";
    break;
    default:
    std::cout << "Not existing (should be between 0 and 2). Set default (1)";
    fResponseType = 1;
  }
  std::cout << std::endl;

  DefineOutput(1, AliMergeableCollection::Class());
}


//________________________________________________________________________
AliAnalysisTaskWeightMTRResponse::~AliAnalysisTaskWeightMTRResponse()
{
  //
  /// Destructor
  //

  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fMergeableCollection;
  }
}

//________________________________________________________________________
void AliAnalysisTaskWeightMTRResponse::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts.SetRun(fInputHandler);
}

//________________________________________________________________________
TH2* AliAnalysisTaskWeightMTRResponse::GetHisto ( const char* source, Bool_t mcResponse, Bool_t perBoard, Bool_t useFit, Int_t makeIt )
{
  /// Get histogram

  TString histoName = Form("weightedHistoPer%s",perBoard?"Board":"Eta");
  if ( useFit ) histoName.Append("_fit");
  TString identifier = Form("/%s/",source);
  identifier += ( mcResponse ) ? "MCRESP" : "DATARESP";
  TH2* histo = static_cast<TH2*>(fMergeableCollection->GetObject(identifier.Data(), histoName.Data()));
  if ( histo || makeIt == 0 ) return histo;

  TString ptTitle = "";

  TString yTitle = perBoard ? "Board ID" : "#eta";
  Int_t nYbins = perBoard ? 234 : 15;
  Double_t yMin = perBoard ? 0.5 : -4.;
  Double_t yMax = perBoard ? 234.5 : -2.5;

  TString mustring = "";
  for ( Int_t imu=0; imu<makeIt; imu++ ) mustring += "#mu";

  histo = new TH2D(histoName,histoName,40,0.,20.,15,-4.,-2.5);
  histo->SetXTitle(Form("p^{%s}_{T} (GeV/c)",mustring.Data()));
  histo->SetYTitle(Form("y_{%s}",mustring.Data()));
  histo->Sumw2();

  fMergeableCollection->Adopt(identifier, histo);
  AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  return histo;
}

//________________________________________________________________________
TString AliAnalysisTaskWeightMTRResponse::GetSingleMuSource ( AliVParticle* track )
{
  /// Get name of single muon source
  if ( fUtilityMuonAncestor.IsUnidentified(track,MCEvent()) ) return "Unidentified";
  if ( fUtilityMuonAncestor.IsHadron(track,MCEvent()) ) return "RecoHadron";
  if ( fUtilityMuonAncestor.IsSecondaryMu(track,MCEvent()) ) return "SecondaryMu";
  if ( fUtilityMuonAncestor.IsDecayMu(track,MCEvent()) ) return "DecayMu";
  if ( fUtilityMuonAncestor.IsBeautyMu(track,MCEvent()) ) return "BeautyMu";
  if ( fUtilityMuonAncestor.IsCharmMu(track,MCEvent()) ) return "CharmMu";
  if ( fUtilityMuonAncestor.IsWBosonMu(track,MCEvent()) ) return "WbosonMu";
  if ( fUtilityMuonAncestor.IsZBosonMu(track,MCEvent()) ) return "ZbosonMu";
  if ( fUtilityMuonAncestor.IsQuarkoniumMu(track,MCEvent()) ) return "QuarkonimuMu";

  return "DecayMu";
}

//___________________________________________________________________________
void AliAnalysisTaskWeightMTRResponse::UserCreateOutputObjects()
{
  //
  /// Create output histograms
  //

  fMergeableCollection = new AliMergeableCollection(GetOutputSlot(1)->GetContainer()->GetName());

  UInt_t mask = fMuonTrackCuts.GetFilterMask();
  mask &= ~(AliMuonTrackCuts::kMuMatchApt | AliMuonTrackCuts::kMuMatchLpt | AliMuonTrackCuts::kMuMatchHpt);
  switch ( fResponseType ) {
    case AliMTRParameterizedResponse::kLptOverApt:
    mask |= AliMuonTrackCuts::kMuMatchApt;
    break;
    case AliMTRParameterizedResponse::kHptOverLpt:
    mask |= AliMuonTrackCuts::kMuMatchLpt;
    break;
  }
  fMuonTrackCuts.SetFilterMask(mask);

  fMuonTrackCuts.Print("mask");

  PostData(1,fMergeableCollection);
}

//________________________________________________________________________
void AliAnalysisTaskWeightMTRResponse::UserExec ( Option_t * /*option*/ )
{

  /// Fill histograms

  if ( ! MCEvent() ) {
    AliError("Need MC!");
    return;
  }

  AliVParticle* track = 0x0;
  AliTrackMore* trackMore = 0x0;

  std::vector<AliTrackMore*> selected;
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(InputEvent());
  for ( Int_t itrack = 0; itrack < nTracks; itrack++ ) {
    track = AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());
    if ( ! fMuonTrackCuts.IsSelected(track) ) continue;
    trackMore = new AliTrackMore(track);
    trackMore->SetParticleType(fUtilityDimuonSource.GetParticleType(track,MCEvent()));
    trackMore->SetHistory(AliAnalysisMuonUtility::GetTrackHistory(track,MCEvent()));
    TString singleMuName = GetSingleMuSource(track);
    Int_t board = ( fResponseType>0 ) ? AliAnalysisMuonUtility::GetLoCircuit(track) : 0;
    Double_t muPt = track->Pt();
    Double_t muEta = track->Eta();
    Double_t muRapidity = track->Y();
    for ( Int_t ifit=0; ifit<2; ifit++ ) {
      for ( Int_t imc=0; imc<2; imc++ ) {
        for ( Int_t iperboard=0; iperboard<2; iperboard++ ) {
          Double_t wgt = ( iperboard == 0 ) ? fResponse.WeightPerEta(muPt,muEta,fResponseType,imc,ifit) : fResponse.WeightPerBoard(muPt,board,fResponseType,imc,ifit);
          trackMore->SetWgt(wgt,iperboard,imc,ifit);
          // printf("Pt %g  Eta %g  Rap. %g  board %i  isMC %i  isFit %i  wgt %g\n",muPt,muEta,muRapidity,board,imc,ifit,wgt); // REMEMBER TO CUT
          if ( wgt != 0. ) static_cast<TH2*>(GetHisto(singleMuName.Data(),imc,iperboard,ifit,1))->Fill(muPt,muRapidity,wgt);
        } // isPerBoard
      } // isMC
    } // isFit
    selected.push_back(trackMore);
  }

  Int_t nSelected = selected.size();
  for ( Int_t imu1=0; imu1<nSelected; imu1++ ) {
    for ( Int_t imu2=imu1+1; imu2<nSelected; imu2++ ) {
      Int_t commonAncestor = fUtilityDimuonSource.GetCommonAncestor(selected[imu1]->GetTrack(),selected[imu2]->GetTrack(),MCEvent());
      TString pairType = fUtilityDimuonSource.GetPairType(selected[imu1]->GetParticleType(), selected[imu2]->GetParticleType(), commonAncestor, MCEvent());

      TLorentzVector dimu = AliAnalysisMuonUtility::GetTrackPair(selected[imu1]->GetTrack(),selected[imu2]->GetTrack());
      Double_t pt = dimu.Pt();
      Double_t rapidity = dimu.Rapidity();
      for ( Int_t ifit=0; ifit<2; ifit++ ) {
        for ( Int_t imc=0; imc<2; imc++ ) {
          for ( Int_t iperboard=0; iperboard<2; iperboard++ ) {
            Double_t wgt = selected[imu1]->GetWgt(iperboard,imc,ifit) * selected[imu2]->GetWgt(iperboard,imc,ifit);
            if ( wgt == 0. ) continue;
            static_cast<TH2*>(GetHisto(pairType.Data(),imc,iperboard,ifit,2))->Fill(pt,rapidity,wgt);
          } // isPerBoard
        } // isMC
      } // isFit
    } // loop on mu2
  } // loop on mu1

  for ( auto& trm : selected ) delete trm;

  PostData(1,fMergeableCollection);
}

//________________________________________________________________________
void AliAnalysisTaskWeightMTRResponse::Terminate(Option_t *)
{
  //
  /// Draw some histograms at the end.
  //

  fMergeableCollection = static_cast<AliMergeableCollection*>(GetOutputData(1));
  if ( ! fMergeableCollection ) return;

  TList* sources = fMergeableCollection->CreateListOfKeys(0);
  TIter next(sources);
  TObject* src = 0x0;
  while ( (src = next()) ) {
    for ( Int_t iperboard=0; iperboard<2; iperboard++ ) {
      for ( Int_t ifit=0; ifit<2; ifit++ ) {
        TH2* histoMC = GetHisto(src->GetName(),kTRUE,iperboard,ifit);
        if ( ! histoMC ) continue;
        TH2* histoData = GetHisto(src->GetName(),kFALSE,iperboard,ifit);
        if ( ! histoData ) continue;
        TString histoName = histoMC->GetName();
        histoName.Prepend(Form("%s_relDiff_",src->GetName()));
        TH1* histoRatio = static_cast<TH1*>(histoMC->Clone(histoName.Data()));
        TString canName = Form("%s_%s",GetName(),histoRatio->GetName());
        TCanvas* can = new TCanvas(canName.Data(),canName.Data(),100*iperboard,0,800,800);
        histoRatio->Add(histoData,-1.);
        histoRatio->Divide(histoData);
        histoName.Append("   (MC/data)-1");
        histoRatio->SetTitle(histoName.Data());
        histoRatio->Draw("colz");
      }
    }
  }
}

////////////////////////////////////////////////////////////////////
//
// AliTrackMore
//
////////////////////////////////////////////////////////////////////

//________________________________________________________________________
AliAnalysisTaskWeightMTRResponse::AliTrackMore::AliTrackMore ( AliVParticle* track ):
TObject(),
fTrack(track),
// fTrigClassCut(0),
fParticleType(-1),
fAncestor(-1),
fHistory("")
{
  /// Ctr
}


//________________________________________________________________________
AliAnalysisTaskWeightMTRResponse::AliTrackMore::~AliTrackMore()
{
  /// Destructor (does nothing since fTrack is not owner)
  fTrack = 0x0;
}

//________________________________________________________________________
void AliAnalysisTaskWeightMTRResponse::AliTrackMore::SetWgt ( Double_t wgt, Bool_t isPerBoard, Bool_t isMC, Bool_t isFit )
{
  /// Set weight
  Int_t iwgt = 4*isPerBoard+2*isFit+isMC;
  fWgts[iwgt] = wgt;
}

//________________________________________________________________________
Double_t AliAnalysisTaskWeightMTRResponse::AliTrackMore::GetWgt ( Bool_t isPerBoard, Bool_t isMC, Bool_t isFit )
{
  /// Get weight
  Int_t iwgt = 4*isPerBoard+2*isFit+isMC;
  return fWgts[iwgt];
}
