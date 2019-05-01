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

#include "AliAnalysisTaskMTRResponse.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMap.h"
#include "THashList.h"
#include "TGraphAsymmErrors.h"

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

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMTRResponse) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskMTRResponse::AliAnalysisTaskMTRResponse() :
  AliAnalysisTaskSE(),
  fMatchTrigKeys(0x0),
  fMergeableCollection(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskMTRResponse::AliAnalysisTaskMTRResponse(const char *name) :
  AliAnalysisTaskSE(name),
  fMatchTrigKeys(0x0),
  fMergeableCollection(0x0)
{
  //
  /// Constructor
  //

  TString matchNames = "MatchNo MatchAllPt MatchLowPt MatchHighPt";
  fMatchTrigKeys = matchNames.Tokenize(" ");

  DefineOutput(1, AliMergeableCollection::Class());
}


//________________________________________________________________________
AliAnalysisTaskMTRResponse::~AliAnalysisTaskMTRResponse()
{
  //
  /// Destructor
  //

  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fMergeableCollection;
  }
  delete fMatchTrigKeys;
}

//________________________________________________________________________
void AliAnalysisTaskMTRResponse::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts.SetRun(fInputHandler);
}

//________________________________________________________________________
TH2* AliAnalysisTaskMTRResponse::GetHisto ( TString identifier, Int_t imatch, Bool_t perBoard, Bool_t makeIt )
{
  /// Get histogram

  TString histoName = Form("histo%s%s",fMatchTrigKeys->At(imatch)->GetName(),perBoard?"PerBoard":"PerEta");

  TH2* histo = static_cast<TH2*>(fMergeableCollection->GetObject(identifier.Data(), histoName.Data()));
  if ( histo || ! makeIt ) return histo;

  TString ptTitle = "p_{T} (GeV/c)";

  // We need large granularity at low pt
  // and decreasing granularity at larger pt
  // We define here the bin width and the bin limits
  Double_t stepRanges[] = {0., 3., 6., 10., 14., 20., 30.};
  const Int_t kNsteps = sizeof(stepRanges)/sizeof(stepRanges[0])-1;
  Double_t steps[kNsteps] = {0.1,0.2,0.5,1.,2.,10.};

  // And we then build the pt axis according to the granularity defined above
  Int_t countBins = 0;
  for ( Int_t istep=0; istep<kNsteps; istep++ ) {
    countBins += TMath::Nint((stepRanges[istep+1]-stepRanges[istep])/steps[istep]);
  }
  const Int_t kNptBins = countBins;
  Double_t ptBins[kNptBins+1];
  ptBins[0] = 0;
  Int_t currStep = 0;
  for ( Int_t ibin=1; ibin<=kNptBins; ibin++ ) {
    if ( ptBins[ibin-1] >= stepRanges[currStep+1]-0.001) currStep++;
    ptBins[ibin] = ptBins[ibin-1]+steps[currStep];
  }

  // The binning in eta cannot be too small
  // otherwise the probed region is smaller than the size of one board
  // a rebin will be needed

  TString yTitle = perBoard ? "Board ID" : "#eta";
  Int_t nYbins = perBoard ? 234 : 15;
  Double_t yMin = perBoard ? 0.5 : -4.;
  Double_t yMax = perBoard ? 234.5 : -2.5;

  histo = new TH2D(histoName,histoName,kNptBins,ptBins,nYbins,yMin,yMax);
  histo->SetXTitle(ptTitle.Data());
  histo->SetYTitle(yTitle.Data());

  fMergeableCollection->Adopt(identifier, histo);
  AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  return histo;
}


//___________________________________________________________________________
void AliAnalysisTaskMTRResponse::FinishTaskOutput()
{
  /// Sum histograms

  TString histoName = "", auxName = "";
  TIter next(fMergeableCollection->Map());

  TObjString* str;

  while ( (str = static_cast<TObjString*>(next())) ) {
    for ( Int_t iresp=0; iresp<2; iresp++ ) {
      Bool_t hasNoMatch = kFALSE;
      for ( Int_t imatch=0; imatch<3; imatch++ ) {
        TH2* histo = GetHisto(str->GetName(),imatch,iresp);
        if ( ! histo ) {
          if ( imatch == 1 && hasNoMatch ) {
            // If Apt is not found, but noMatch is, it is because Lpt=Apt=0.5 GeV/c
            // In this case, fill the Apt case as well
            histo = GetHisto(str->GetName(),imatch,iresp,kTRUE);
          }
          else continue;
        }
        if ( imatch == 0 ) hasNoMatch = kTRUE;
        for ( Int_t jmatch=imatch+1; jmatch<4; jmatch++ ) {
          TH2* auxHisto = GetHisto(str->GetName(),jmatch,iresp);
          if ( ! auxHisto ) continue;
          histo->Add(auxHisto);
        } // loop on higher pt match histograms
      } // loop on histograms
    }
  } // loop on identifiers
//    THashList* list = static_cast<THashList*>(fMergeableCollection->Map()->GetValue(str->GetName()));
//    TIter nextHisto(list);
//    while ( (histo = static_cast<TH2*>(nextHisto())) ) {
//      histoName = histo->GetName();
//      Int_t trigMatch = -1;
//      for ( Int_t imatch=0; imatch<3; imatch++ ) {
//        if ( histoName.Contains(fMatchTrigKeys->UncheckedAt(imatch)->GetName()) ) {
//          trigMatch = imatch;
//          break;
//        }
//      }
//
//      if ( trigMatch < 0 ) continue;
//
//      for ( Int_t imatch=trigMatch+1; imatch<4; imatch++ ) {
//        auxName = histoName;
//        auxName.ReplaceAll(fMatchTrigKeys->UncheckedAt(trigMatch)->GetName(),fMatchTrigKeys->UncheckedAt(imatch)->GetName());
//        TH2* auxHisto = static_cast<TH2*>(list->FindObject(auxName.Data()));
//        if ( ! auxHisto ) continue;
//        histo->Add(auxHisto);
//      }
//    } // loop on histograms
//  } // loop on identifiers
}

//___________________________________________________________________________
void AliAnalysisTaskMTRResponse::SetMuonEventCuts ( AliMuonEventCuts* muonEventCuts, const char* mbClassPattern, const char* muLptClassPattern )
{
  /// Set muon event cuts for muon triggered and MB events
  fMuonEventCutsMuTrig = *muonEventCuts;
  fMuonEventCutsMB = *muonEventCuts;

  fMuonEventCutsMuTrig.SetTrigClassPatterns(muLptClassPattern,"");
  fMuonEventCutsMB.SetTrigClassPatterns(mbClassPattern,"");
}

//___________________________________________________________________________
void AliAnalysisTaskMTRResponse::UserCreateOutputObjects()
{
  //
  /// Create output histograms
  //

  fMergeableCollection = new AliMergeableCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fMuonEventCutsMB.Print("mask");
  fMuonTrackCuts.Print("mask");

  PostData(1,fMergeableCollection);
}

//________________________________________________________________________
void AliAnalysisTaskMTRResponse::UserExec ( Option_t * /*option*/ )
{

  /// Fill histograms

  Bool_t isMC = kFALSE;
  if ( MCEvent() ) isMC = kTRUE;

  Bool_t matchTrigClass[2];
  matchTrigClass[0] = fMuonEventCutsMB.IsSelected(fInputHandler);
  matchTrigClass[1] = isMC ? kFALSE : fMuonEventCutsMuTrig.IsSelected(fInputHandler);
  if ( ! matchTrigClass[0] && ! matchTrigClass[1] ) return;

  AliVParticle* track = 0x0;
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(InputEvent());
  for ( Int_t itrack = 0; itrack < nTracks; itrack++ ) {
    track = AliAnalysisMuonUtility::GetTrack(itrack,InputEvent());
    if ( ! fMuonTrackCuts.IsSelected(track) ) continue;

    if ( isMC ) {
      // Select muons
      if ( track->GetLabel() < 0 ) continue;
      AliVParticle* mcTrack = MCEvent()->GetTrack(track->GetLabel());
      if ( TMath::Abs(mcTrack->PdgCode()) != 13 ) continue;
    }

    Int_t matchTrig = AliAnalysisMuonUtility::GetMatchTrigger(track);
    Double_t pt = track->Pt();
    Int_t board = AliAnalysisMuonUtility::GetLoCircuit(track);

    for ( Int_t itrig=0; itrig<2; itrig++ ) {
      if ( ! matchTrigClass[itrig] ) continue;
      if ( itrig == 1 && matchTrig < 2 ) continue;
      TString identifier = ( itrig == 0 ) ? "MB" : "MuTrig";
      if ( board > 0 ) GetHisto(identifier,matchTrig,kTRUE,kTRUE)->Fill(pt,board);
      GetHisto(identifier,matchTrig,kFALSE,kTRUE)->Fill(pt,track->Eta());
    }
  }

  PostData(1,fMergeableCollection);
}


//________________________________________________________________________
void AliAnalysisTaskMTRResponse::Terminate(Option_t *)
{
  //
  /// Draw some histograms at the end.
  //

  fMergeableCollection = static_cast<AliMergeableCollection*>(GetOutputData(1));
  if ( ! fMergeableCollection ) return;

  Bool_t isMC = kFALSE;

  TH2* histos[4];
  for ( Int_t itrig=0; itrig<2; itrig++ ) {
    TString identifier = ( itrig == 0 ) ? "MB" : "MuTrig";
    for ( Int_t imatch=0; imatch<4; imatch++ ) {
      histos[imatch] = GetHisto(identifier,imatch,kFALSE,kFALSE);
    }

    Int_t firstMatch = 2;
    Int_t lastMatch = 3;
    if ( histos[0] ) {
      // This is MC
      firstMatch = 1;
      lastMatch = 2;
      isMC = kTRUE;
    }
    else if ( itrig == 1 ) {
      firstMatch = 3;
      lastMatch = 3;
    }

    for ( Int_t imatch=firstMatch; imatch<=lastMatch; imatch++ ) {
      if ( ! histos[imatch] || ! histos[imatch-1] ) continue;
      TString canName = Form("can%s%sOver%s",GetName(),fMatchTrigKeys->At(imatch)->GetName(),fMatchTrigKeys->At(imatch-1)->GetName());
      TCanvas* can = new TCanvas(canName.Data(),canName.Data(),0,0,1200,800);
      TAxis* axis = histos[imatch]->GetYaxis();
      if ( axis->GetNbins() == 15 ) can->Divide(5,3);
      else can->DivideSquare(axis->GetNbins());
      for ( Int_t ibin=1; ibin<=axis->GetNbins(); ibin++ ) {
        TH1* projNum = histos[imatch]->ProjectionX("tmpNum",ibin,ibin);
        TH1* projDen = histos[imatch-1]->ProjectionX("tmpDen",ibin,ibin);
        TGraphAsymmErrors* gr = new TGraphAsymmErrors(projNum,projDen);
        gr->SetTitle(Form("%g<%s<%g",axis->GetBinLowEdge(ibin),axis->GetTitle(),axis->GetBinUpEdge(ibin)));
        gr->GetXaxis()->SetTitle(histos[imatch]->GetXaxis()->GetTitle());
        gr->GetXaxis()->SetRangeUser(0.,imatch==3?12.:4.);
        delete projNum;
        delete projDen;
        can->cd(ibin);
        gr->Draw("ap");
      }
    } // loop on numerator
    if ( isMC ) break;
  } // loop on trigger classes
}
