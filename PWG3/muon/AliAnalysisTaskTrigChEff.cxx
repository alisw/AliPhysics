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

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskSingleMu
/// Analysis task for single muons in the spectrometer.
/// The output is a list of histograms.
/// The macro class can run on AOD or ESDs.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//    Implementation of the class for trigger chamber efficiency determinaltion
//----------------------------------------------------------------------------


#define AliAnalysisTaskTrigChEff_cxx

// ROOT includes
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"

// STEER includes
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

// ANALYSIS includes
#include "AliAnalysisTaskSE.h"

#include "AliAnalysisTaskTrigChEff.h"

ClassImp(AliAnalysisTaskTrigChEff)

//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff(const char *name) :
  AliAnalysisTaskSE(name), 
  fUseGhosts(kFALSE),
  fList(0)
{
  //
  /// Constructor.
  //
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,  TList::Class());
}

//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::UserCreateOutputObjects() {
  //
  /// Create histograms
  /// Called once
  //

  TString cathCode[2] = {"bendPlane", "nonBendPlane"};
  TString countTypeName[2] = {"CountInCh", "NonCountInCh"};

  const Char_t* yAxisTitle = "counts";

  const Int_t kNboards = 234; //AliMpConstants::NofLocalBoards();
  const Int_t kFirstTrigCh = 11;//AliMpConstants::NofTrackingChambers()+1;

  Int_t chamberBins = kNchambers;
  Float_t chamberLow = kFirstTrigCh-0.5, chamberHigh = kFirstTrigCh+kNchambers-0.5;
  const Char_t* chamberName = "chamber";

  Int_t slatBins = kNslats;
  Float_t slatLow = 0-0.5, slatHigh = kNslats-0.5;
  const Char_t* slatName = "slat";

  Int_t boardBins = kNboards;
  Float_t boardLow = 1-0.5, boardHigh = kNboards+1.-0.5;
  const Char_t* boardName = "board";

  Int_t angleBins = 280;
  Float_t angleLow = -70., angleHigh = 70.;
  const Char_t* angleNameX = "#theta_{x} (deg)";
  const Char_t* angleNameY = "#theta_{y} (deg)";

  TString baseName, histoName;
  fList = new TList();

  TH1F* histo;

  histo = new TH1F("nTracksInSlat", "Num. of tracks used for efficiency calculation", 
		   slatBins, slatLow, slatHigh);
  histo->GetXaxis()->SetTitle(slatName);
  histo->GetYaxis()->SetTitle("num of used tracks");

  fList->AddAt(histo, kHtracksInSlat);

  histo = new TH1F("nTracksInBoard", "Num. of tracks used for efficiency calculation", 
		   boardBins, boardLow, boardHigh);
  histo->GetXaxis()->SetTitle(boardName);
  histo->GetYaxis()->SetTitle("num of used tracks");

  fList->AddAt(histo, kHtracksInBoard);

  for(Int_t hType=0; hType<kNcounts; hType++){
    Int_t hindex = (hType==0) ? kHchamberEff : kHchamberNonEff;
    for(Int_t cath=0; cath<kNcathodes; cath++){
      histoName = Form("%sChamber%s", cathCode[cath].Data(), countTypeName[hType].Data());
      histo = new TH1F(histoName, histoName,
		       chamberBins, chamberLow, chamberHigh);
      histo->GetXaxis()->SetTitle(chamberName);
      histo->GetYaxis()->SetTitle(yAxisTitle);
	
      fList->AddAt(histo, hindex + cath);
    } // loop on cath
  } // loop on counts

  for(Int_t hType=0; hType<kNcounts; hType++){
    Int_t hindex = (hType==0) ? kHslatEff : kHslatNonEff;
    for(Int_t cath=0; cath<kNcathodes; cath++){
      for(Int_t ch=0; ch<kNchambers; ch++){
	Int_t chCath = GetPlane(cath, ch);
	histoName = Form("%sSlat%s%i", cathCode[cath].Data(), countTypeName[hType].Data(), kFirstTrigCh+ch);
	histo = new TH1F(histoName, histoName,
			 slatBins, slatLow, slatHigh);
	histo->GetXaxis()->SetTitle(slatName);
	histo->GetYaxis()->SetTitle(yAxisTitle);

	fList->AddAt(histo, hindex + chCath);
      } // loop on chamber
    } // loop on cath
  } // loop on counts

  for(Int_t hType=0; hType<kNcounts; hType++){
    Int_t hindex = (hType==0) ? kHboardEff : kHboardNonEff;
    for(Int_t cath=0; cath<kNcathodes; cath++){
      for(Int_t ch=0; ch<kNchambers; ch++){
	Int_t chCath = GetPlane(cath, ch);
	histoName = Form("%sBoard%s%i", cathCode[cath].Data(), countTypeName[hType].Data(), kFirstTrigCh+ch);
	histo = new TH1F(histoName, histoName,
			 boardBins, boardLow, boardHigh);
	histo->GetXaxis()->SetTitle(boardName);
	histo->GetYaxis()->SetTitle(yAxisTitle);

	fList->AddAt(histo, hindex + chCath);
      } // loop on chamber
    } // loop on cath
  } // loop on counts

  histo = new TH1F("thetaX", "Angular distribution",
		   angleBins, angleLow, angleHigh);
  histo->GetXaxis()->SetTitle(angleNameX);
  histo->GetYaxis()->SetTitle("entries");
  fList->AddAt(histo, kHthetaX);

  histo = new TH1F("thetaY", "Angular distribution",
		   angleBins, angleLow, angleHigh);
  histo->GetXaxis()->SetTitle(angleNameY);
  histo->GetYaxis()->SetTitle("entries");
  fList->AddAt(histo, kHthetaY);

}

//________________________________________________________________________
void AliAnalysisTaskTrigChEff::UserExec(Option_t *) {
  //
  /// Main loop
  /// Called for each event
  //
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*> (InputEvent());

  if (!esdEvent) {
    Printf("ERROR: esdEvent not available\n");
    return;
  }

  Int_t slat = 0, board = 0;
  UShort_t pattern = 0;
  AliESDMuonTrack *esdTrack = 0x0;

  const Float_t kRadToDeg = 180./TMath::Pi();
  Int_t nTracks = esdEvent->GetNumberOfMuonTracks();

  const Int_t kFirstTrigCh = 11; //AliMpConstants::NofTrackingChambers()+1;

  TArrayI othersEfficient(kNchambers);

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    esdTrack = esdEvent->GetMuonTrack(itrack);

    if ( ! esdTrack->ContainTrackerData() && ! fUseGhosts ) continue;

    pattern =  esdTrack->GetHitsPatternInTrigCh();
    Int_t effFlag = AliESDMuonTrack::GetEffFlag(pattern);

    if(effFlag < AliESDMuonTrack::kChEff) continue; // Track not good for efficiency calculation

    ((TH1F*)fList->At(kHthetaX))->Fill(esdTrack->GetThetaX() * kRadToDeg);
    ((TH1F*)fList->At(kHthetaY))->Fill(esdTrack->GetThetaY() * kRadToDeg);

    othersEfficient.Reset(1);
    for(Int_t cath=0; cath<kNcathodes; cath++){
      for(Int_t ich=0; ich<kNchambers; ich++){
	if( ! AliESDMuonTrack::IsChamberHit(pattern, cath, ich)){
	  for(Int_t jch=0; jch<kNchambers; jch++){
	    if ( jch != ich) {
	      othersEfficient[jch] = 0;
	      //AliInfo(Form("%s ch %i by New", baseOutString.Data(), jch));
	    }
	  } // loop on other chambers
	  break;
	} // if chamber not efficient
      } // loop on chambers
    } // loop on cathodes

    Bool_t rejectTrack = kTRUE;
    for (Int_t ich=0; ich<kNchambers; ich++){
      if ( othersEfficient[ich] > 0 ){
	rejectTrack = kFALSE;
	break;
      }
    }

    if ( rejectTrack ) continue;

    slat = AliESDMuonTrack::GetSlatOrInfo(pattern);
    board = esdTrack->LoCircuit();

    if(effFlag >= AliESDMuonTrack::kSlatEff) ((TH1F*)fList->At(kHtracksInSlat))->Fill(slat);
    if(effFlag >= AliESDMuonTrack::kBoardEff) ((TH1F*)fList->At(kHtracksInBoard))->Fill(board);

    for(Int_t cath=0; cath<kNcathodes; cath++){
      for(Int_t ch=0; ch<kNchambers; ch++){
	if ( ! othersEfficient[ch] )
	  continue; // Reject track if the info of the chamber under study 
	            // is necessary to create the track itself

	Int_t whichType = AliESDMuonTrack::IsChamberHit(pattern, cath, ch) ? kChHit : kChNonHit;

	Int_t iChamber = kFirstTrigCh + ch;
	Int_t hindex = ( whichType == kChHit ) ? kHchamberEff : kHchamberNonEff;
	((TH1F*)fList->At(hindex + cath))->Fill(iChamber);

	if(effFlag < AliESDMuonTrack::kSlatEff) continue; // Track crossed different slats
	Int_t chCath = GetPlane(cath, ch);
	hindex = ( whichType == kChHit ) ? kHslatEff : kHslatNonEff;
	((TH1F*)fList->At(hindex + chCath))->Fill(slat);

	if(effFlag < AliESDMuonTrack::kBoardEff) continue; // Track crossed different boards
	hindex = ( whichType == kChHit ) ? kHboardEff : kHboardNonEff;
	((TH1F*)fList->At(hindex + chCath))->Fill(board);
      } // loop on chambers
    } // loop on cathodes
  } // loop on tracks

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fList);
}

//________________________________________________________________________
void AliAnalysisTaskTrigChEff::Terminate(Option_t *) {
  //
  /// Draw result to the screen
  /// Called once at the end of the query.
  //
  if (!gROOT->IsBatch()) {
    TCanvas *can[kNcathodes];
    TH1F *num = 0x0;
    TH1F *den = 0x0;
    for(Int_t cath=0; cath<kNcathodes; cath++){
      TString canName = Form("can%i",cath);
      can[cath] = new TCanvas(canName.Data(),canName.Data(),10*(1+cath),10*(1+cath),310,310);
      can[cath]->SetFillColor(10); can[cath]->SetHighLightColor(10);
      can[cath]->SetLeftMargin(0.15); can[cath]->SetBottomMargin(0.15);  
      can[cath]->Divide(2,2);
      for(Int_t ch=0; ch<kNchambers; ch++){
	Int_t chCath = GetPlane(cath, ch);
	num = (TH1F*)(fList->At(kHboardEff + chCath)->Clone());
	den = (TH1F*)(fList->At(kHboardNonEff + chCath)->Clone());
	den->Add(num);
	num->Divide(den);
	can[cath]->cd(ch+1);
	num->DrawCopy("E");
      }
    }
  }
}
