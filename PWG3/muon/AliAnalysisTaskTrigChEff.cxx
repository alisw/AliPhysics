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
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraphAsymmErrors.h"

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

//________________________________________________________________________
AliAnalysisTaskTrigChEff::~AliAnalysisTaskTrigChEff()
{
  delete fList;
}


//_________________________________________________________________________
void AliAnalysisTaskTrigChEff::UserCreateOutputObjects() {
  //
  /// Create histograms
  /// Called once
  //

  TString countTypeName[kNcounts] = {"bendPlane", "nonBendPlane","bothPlanes", "allTracks"};

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

  TString baseName, histoName, histoTitle;
  fList = new TList();

  TH1F* histo;
  TH2F* histo2D;

  Int_t histoIndex = -1;

  for(Int_t icount=0; icount<kNcounts; icount++){
    histoName = Form("%sCountChamber", countTypeName[icount].Data());
    histo = new TH1F(histoName, histoName,
		     chamberBins, chamberLow, chamberHigh);
    histo->GetXaxis()->SetTitle(chamberName);
    histo->GetYaxis()->SetTitle(yAxisTitle);
    histoIndex = GetHistoIndex(kHchamberEff, icount);
    fList->AddAt(histo, histoIndex);
  } // loop on counts

  for(Int_t icount=0; icount<kNcounts; icount++){
    for(Int_t ch=0; ch<kNchambers; ch++){
      histoName = Form("%sCountSlatCh%i", countTypeName[icount].Data(), kFirstTrigCh+ch);
      histo = new TH1F(histoName, histoName,
		       slatBins, slatLow, slatHigh);
      histo->GetXaxis()->SetTitle(slatName);
      histo->GetYaxis()->SetTitle(yAxisTitle);
      histoIndex = GetHistoIndex(kHslatEff, icount, ch);
      fList->AddAt(histo, histoIndex);
    } // loop on chamber
  } // loop on counts

  for(Int_t icount=0; icount<kNcounts; icount++){
    for(Int_t ch=0; ch<kNchambers; ch++){
      histoName = Form("%sCountBoardCh%i", countTypeName[icount].Data(), kFirstTrigCh+ch);
      histo = new TH1F(histoName, histoName,
		       boardBins, boardLow, boardHigh);
      histo->GetXaxis()->SetTitle(boardName);
      histo->GetYaxis()->SetTitle(yAxisTitle);
      histoIndex = GetHistoIndex(kHboardEff, icount, ch);
      fList->AddAt(histo, histoIndex);
    } // loop on chamber
  } // loop on counts

  histo2D = new TH2F("checkRejectedBoard", "Rejected tracks motivation", 
		     4, 20.5, 24.5, boardBins, boardLow, boardHigh);
  histo2D->GetXaxis()->SetBinLabel(1,"Many pads");
  histo2D->GetXaxis()->SetBinLabel(2,"Few pads");
  histo2D->GetXaxis()->SetBinLabel(3,"Outside geom");
  histo2D->GetXaxis()->SetBinLabel(4,"Tracker track");
  histo2D->GetYaxis()->SetTitle(boardName);
  histoIndex = GetHistoIndex(kHcheckBoard);
  fList->AddAt(histo2D, histoIndex);
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

  Int_t nTracks = esdEvent->GetNumberOfMuonTracks();

  const Int_t kFirstTrigCh = 11; //AliMpConstants::NofTrackingChambers()+1;

  TArrayI othersEfficient(kNchambers);
  Int_t histoIndex = -1;

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    esdTrack = esdEvent->GetMuonTrack(itrack);

    if ( ! esdTrack->ContainTrackerData() && ! fUseGhosts ) continue;

    pattern =  esdTrack->GetHitsPatternInTrigCh();
    Int_t effFlag = AliESDMuonTrack::GetEffFlag(pattern);

    board = esdTrack->LoCircuit();

    if(effFlag < AliESDMuonTrack::kChEff) {
      histoIndex = GetHistoIndex(kHcheckBoard);
      ((TH2F*)fList->At(histoIndex))->Fill(AliESDMuonTrack::GetSlatOrInfo(pattern), board);
      continue; // Track not good for efficiency calculation
    }

    othersEfficient.Reset(1);
    for(Int_t cath=0; cath<kNcathodes; cath++){
      for(Int_t ich=0; ich<kNchambers; ich++){
	if( ! AliESDMuonTrack::IsChamberHit(pattern, cath, ich)){
	  for(Int_t jch=0; jch<kNchambers; jch++){
	    if ( jch != ich) {
	      othersEfficient[jch] = 0;
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

    for(Int_t ch=0; ch<kNchambers; ch++){
      if ( ! othersEfficient[ch] )
	continue; // Reject track if the info of the chamber under study 
                  // is necessary to create the track itself

      Int_t iChamber = kFirstTrigCh + ch;

      Bool_t hitsBend = AliESDMuonTrack::IsChamberHit(pattern, 0, ch);
      Bool_t hitsNonBend = AliESDMuonTrack::IsChamberHit(pattern, 1, ch);

      Bool_t fillHisto[kNcounts] = {
	hitsBend,
	hitsNonBend,
	( hitsBend && hitsNonBend ),
	kTRUE
      };

      for (Int_t icount=0; icount<kNcounts; icount++){
	if ( ! fillHisto[icount] ) continue;

	histoIndex = GetHistoIndex(kHchamberEff, icount);
	((TH1F*)fList->At(histoIndex))->Fill(iChamber);

	if(effFlag < AliESDMuonTrack::kSlatEff) continue; // Track crossed different slats
	histoIndex = GetHistoIndex(kHslatEff, icount, ch);
	((TH1F*)fList->At(histoIndex))->Fill(slat);

	if(effFlag < AliESDMuonTrack::kBoardEff) continue; // Track crossed different boards
	histoIndex = GetHistoIndex(kHboardEff, icount, ch);
	((TH1F*)fList->At(histoIndex))->Fill(board);
      } // loop on chambers
    } // loop on count types
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
  if ( gROOT->IsBatch() ) return;
  fList = dynamic_cast<TList*> (GetOutputData(1));

  if (!fList) return;

  TCanvas *can;
  TH1F *num = 0x0;
  TH1F *den = 0x0;
  TGraphAsymmErrors* effGraph = 0x0;
  TString baseName[3] = {"Chamber", "RPC", "Board"};
  Int_t baseIndex[3] = {kHchamberEff, kHslatEff, kHboardEff};
  TString effName[kNcounts-1] = {"BendPlane", "NonBendPlane", "BothPlanes"};
  Int_t histoIndexNum = -1, histoIndexDen = -1;
  for (Int_t itype=0; itype<3; itype++) {
    for(Int_t icount=0; icount<kNcounts; icount++){
      TString canName = Form("efficiencyPer%s_%s",baseName[itype].Data(),effName[icount].Data());
      can = new TCanvas(canName.Data(),canName.Data(),10*(1+kNcounts*itype+icount),10*(1+kNcounts*itype+icount),310,310);
      can->SetFillColor(10); can->SetHighLightColor(10);
      can->SetLeftMargin(0.15); can->SetBottomMargin(0.15);  
      if ( itype > 0 )
	can->Divide(2,2);

      for(Int_t ch=0; ch<kNchambers; ch++){
	histoIndexNum = GetHistoIndex(baseIndex[itype], icount, ch);
	histoIndexDen = GetHistoIndex(baseIndex[itype], kAllTracks, ch);
	num = (TH1F*)(fList->At(histoIndexNum));
	den = (TH1F*)(fList->At(histoIndexDen));
	effGraph = new TGraphAsymmErrors(num, den);
	effGraph->GetYaxis()->SetRangeUser(0., 1.1);
	effGraph->GetYaxis()->SetTitle("Efficiency");
	effGraph->GetXaxis()->SetTitle(baseName[itype].Data());
	can->cd(ch+1);
	effGraph->Draw("AP");
	if ( itype == 0 ) break;
      } // loop on chamber
    } // loop on count types
  } // loop on histo
}

//________________________________________________________________________
Int_t
AliAnalysisTaskTrigChEff::GetHistoIndex(Int_t histoType, Int_t countType, 
					Int_t chamber)
{
  //
  /// Return the index of the histogram in the list
  //
  switch ( histoType ) {
  case kHchamberEff:
    return 0 + countType;
  case kHslatEff:
    return 4 + kNchambers*countType + chamber;
  case kHboardEff:
    return 20 + kNchambers*countType + chamber;
  case kHcheckBoard:
    return 36;
  }

  /*
  const Int_t kNhistosPlaneCorr = 38;

  switch ( histoType ){
  case kHtracksInSlat:
    return 0 + planeCorrelation*kNhistosPlaneCorr;
  case kHtracksInBoard:
    return 1 + planeCorrelation*kNhistosPlaneCorr;
  case kHchamberEff:
    return 2 + kNcathodes*countType + cathode
      + planeCorrelation*kNhistosPlaneCorr;
  case kHslatEff:
    return 6 + kNchambers*kNcathodes*countType 
      + kNchambers*cathode + chamber
      + planeCorrelation*kNhistosPlaneCorr;
  case kHboardEff:
    return 22 + kNchambers*kNcathodes*countType 
      + kNchambers*cathode + chamber
      + planeCorrelation*kNhistosPlaneCorr;
  case kHcheckBoard:
    return 0 + 2*kNhistosPlaneCorr;
  }
  */
  return -1;
}
