#define AliAnalysisTaskTrigChEff_cxx

// ROOT includes
#include "TChain.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"

// STEER includes
#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"

// ANALYSIS includes
#include "AliAnalysisTask.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskTrigChEff.h"

ClassImp(AliAnalysisTaskTrigChEff)

//________________________________________________________________________
AliAnalysisTaskTrigChEff::AliAnalysisTaskTrigChEff(const char *name) :
  AliAnalysisTask(name,""), 
  fESD(0),
  fAOD(0),
  fAnalysisType("ESD"),
  fList(0)
{
  //
  /// Constructor.
  //
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TObjArray container
  DefineOutput(0,  TList::Class());
}

//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::ConnectInputData(Option_t *) {
  //
  /// Connect ESD or AOD here
  /// Called once
  //

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    if(fAnalysisType == "ESD") {
      tree->SetBranchStatus("*", kFALSE);
      tree->SetBranchStatus("MuonTracks.*", kTRUE);

      AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
      if (!esdH) {
	Printf("ERROR: Could not get ESDInputHandler");
      } else
	fESD = esdH->GetEvent();
    }
    else if(fAnalysisType == "AOD") {
      AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      
      if (!aodH) {
	Printf("ERROR: Could not get AODInputHandler");
      } else
	fAOD = aodH->GetEvent();
    }
    else 
      Printf("Wrong analysis type: Only ESD and AOD types are allowed!");
  }
}

//___________________________________________________________________________
void AliAnalysisTaskTrigChEff::CreateOutputObjects() {
  //
  /// Create histograms
  /// Called once
  //
  Printf("   CreateOutputObjects of task %s\n", GetName());

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
    Int_t hindex = (hType==0) ? kHchamberAllEff : kHchamberNonEff;
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
    Int_t hindex = (hType==0) ? kHslatAllEff : kHslatNonEff;
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
    Int_t hindex = (hType==0) ? kHboardAllEff : kHboardNonEff;
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
}

//________________________________________________________________________
void AliAnalysisTaskTrigChEff::Exec(Option_t *) {
  //
  /// Main loop
  /// Called for each event
  //
  Int_t nTracks = 0, board = 0;
  UShort_t pattern = 0;
  AliESDMuonTrack *esdTrack = 0x0;
  AliAODTrack* aodTrack = 0x0;

  if(fAnalysisType == "ESD") {
    if (!fESD) {
      Printf("ERROR: fESD not available");
      return;
    }
    nTracks = fESD->GetNumberOfMuonTracks(); 
  }
  else if(fAnalysisType == "AOD") {
    if (!fAOD) {
      Printf("ERROR: fAOD not available");
      return;
    }
    nTracks = fAOD->GetNumberOfTracks();
  }

  // Object declaration
  const Int_t kFirstTrigCh = 11; //AliMpConstants::NofTrackingChambers()+1;

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {
    if(fAnalysisType == "ESD") {
      esdTrack = fESD->GetMuonTrack(itrack);
      pattern =  esdTrack->GetHitsPatternInTrigCh();
      board = esdTrack->LoCircuit();
    }
    else if(fAnalysisType == "AOD") {
      aodTrack = fAOD->GetTrack(itrack);
      if(!aodTrack->IsMuonTrack()) continue;
      pattern =  aodTrack->GetHitsPatternInTrigCh();
      board = 0; // aodTrack->LoCircuit(); Lo Circuit not implemented in AOD
    }

    Int_t effFlag = GetEffFlag(pattern);

    if(effFlag < kChEff) continue; // Track not good for efficiency calculation

    Int_t slat = GetSlat(pattern);

    if(effFlag >= kSlatEff) ((TH1F*)fList->At(kHtracksInSlat))->Fill(slat);
    if(effFlag >= kBoardEff) ((TH1F*)fList->At(kHtracksInBoard))->Fill(board);

    for(Int_t cath=0; cath<kNcathodes; cath++){
      Int_t ineffCh = IsChInefficient(pattern, cath);
      Int_t nChambers = kNchambers;
      for(Int_t ch=0; ch<nChambers; ch++){
	Int_t whichType = kAllChEff;
	Int_t currCh = ch;
	if(ineffCh>=0){
	  whichType = kChNonEff;
	  currCh = ineffCh;
	  nChambers = -1;
	}

	Int_t iChamber = kFirstTrigCh + currCh;
	Int_t hindex = (whichType==kAllChEff) ? kHchamberAllEff : kHchamberNonEff;
	((TH1F*)fList->At(hindex + cath))->Fill(iChamber);

	if(effFlag < kSlatEff) continue; // Track crossed different slats
	Int_t chCath = GetPlane(cath, currCh);
	hindex = (whichType==kAllChEff) ? kHslatAllEff : kHslatNonEff;
	((TH1F*)fList->At(hindex + chCath))->Fill(slat);

	if(effFlag < kBoardEff) continue; // Track crossed different boards
	hindex = (whichType==kAllChEff) ? kHboardAllEff : kHboardNonEff;
	((TH1F*)fList->At(hindex + chCath))->Fill(board);
      } // loop on chambers
    } // loop on cathodes
  }

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fList);
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
	num = (TH1F*)(fList->At(kHboardAllEff + chCath)->Clone());
	den = (TH1F*)(fList->At(kHboardNonEff + chCath)->Clone());
	den->Add(num);
	num->Divide(den);
	can[cath]->cd(ch+1);
	num->DrawCopy("E");
      }
    }
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskTrigChEff::IsChInefficient(UShort_t pattern, 
						Int_t cathode)
{
  //
  /// Check which chamber was inefficient.
  //
  Int_t ineffCh = -999;
  for(Int_t ch=0; ch<kNchambers; ch++){
    Int_t chCath = GetPlane(cathode, ch);
    Int_t invert = kNplanes - chCath - 1;
    Int_t response = (pattern >> invert) & 0x01;
    if(!response) ineffCh = ch;
  }
  return ineffCh;
}
