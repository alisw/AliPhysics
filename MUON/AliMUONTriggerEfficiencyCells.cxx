/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// $Id$

#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMpConstants.h"

// Classes for display
#include "AliMUONTriggerDisplay.h"
#include "AliCDBManager.h"
#include "AliMpDDLStore.h"
#include "AliMpLocalBoard.h"
#include "AliMpPad.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

#include "AliLog.h"

#include "TRandom.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TCanvas.h"

#include <fstream>
#include <cassert>

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerEfficiencyCells
/// A class to store and give access to the trigger chamber efficiency.
///
/// Efficiency is stored per cathode on local boards
///
/// The main method of this class is IsTriggered().
///
/// $ALICE_ROOT/MUON/data/efficiencyCells.dat contains efficiency 
/// for each chamber (i.e. DetElement). 
///
/// In the case of local boards, efficiency is stored from left to right
/// per increasing board number (from 1 to 234)
///
/// The file can be edited in order to change efficiency
/// in a chosen local board/region of the chamber.
///
///
/// But please note that this object is also available from the CDB 
///
/// \author Diego Stocco; INFN Torino
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONTriggerEfficiencyCells)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells()
:
TObject(),
fCountHistoList(0x0),
fNoCountHistoList(0x0),
fFiredStrips(0x0),
fDisplayHistoList(0x0),
fBoardLabelList(0x0),
fFiredFitHistoList(0x0),
fFiredDisplayHistoList(0x0)
{
///  Default constructor.
  CheckConstants();
  Reset();
  InitHistos();
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const Char_t* filename)
:
TObject(),
fCountHistoList(0x0),
fNoCountHistoList(0x0),
fFiredStrips(0x0),
fDisplayHistoList(0x0),
fBoardLabelList(0x0),
fFiredFitHistoList(0x0),
fFiredDisplayHistoList(0x0)
{
///  Constructor using an ASCII file.
  CheckConstants();
  Reset();
  ReadFile(filename);
}

AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(TList *countHistoList,
							     TList *noCountHistoList)
:
TObject(),
fCountHistoList(countHistoList),
fNoCountHistoList(noCountHistoList),
fFiredStrips(0x0),
fDisplayHistoList(0x0),
fBoardLabelList(0x0),
fFiredFitHistoList(0x0),
fFiredDisplayHistoList(0x0)
{
///  Constructor using an ASCII file.
  CheckConstants();
  Reset();
  InitHistos();
  FillHistosFromList();
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const AliMUONTriggerEfficiencyCells& other)
:
TObject(other),
fCountHistoList(other.fCountHistoList),
fNoCountHistoList(other.fNoCountHistoList),
fFiredStrips(other.fFiredStrips),
fDisplayHistoList(other.fDisplayHistoList),
fBoardLabelList(other.fBoardLabelList),
fFiredFitHistoList(other.fFiredFitHistoList),
fFiredDisplayHistoList(other.fFiredDisplayHistoList)
{
/// Copy constructor

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    fBoardEfficiency[chCath] = other.fBoardEfficiency[chCath];
    fSlatEfficiency[chCath] = other.fSlatEfficiency[chCath];
  }
}

//_____________________________________________________________________________
AliMUONTriggerEfficiencyCells& AliMUONTriggerEfficiencyCells::operator=(const AliMUONTriggerEfficiencyCells& other)
{
  /// Asignment operator
  // check assignement to self
  if (this == &other)
    return *this;

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    fBoardEfficiency[chCath] = other.fBoardEfficiency[chCath];
    fSlatEfficiency[chCath] = other.fSlatEfficiency[chCath];
  }

  fCountHistoList = other.fCountHistoList;
  fNoCountHistoList = other.fNoCountHistoList;
  fFiredStrips = other.fFiredStrips;

  fDisplayHistoList = other.fDisplayHistoList;
  fBoardLabelList = other.fBoardLabelList;
  fFiredFitHistoList = other.fFiredFitHistoList;
  fFiredDisplayHistoList = other.fFiredDisplayHistoList;
    
  return *this;
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::~AliMUONTriggerEfficiencyCells()
{
///  Destructor.
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::GetCellEfficiency(Int_t detElemId, Int_t localBoard, Float_t &eff1, Float_t &eff2) const
{
///  Get the efficiencies of the 2 cathodes at a given local board

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t bin = fBoardEfficiency[chamber]->FindBin(localBoard);
  eff1 = fBoardEfficiency[chamber]->GetBinContent(bin);
  eff2 = fBoardEfficiency[fgkNchambers+chamber]->GetBinContent(bin);
}


//__________________________________________________________________________
void 
AliMUONTriggerEfficiencyCells::IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trig1, Bool_t &trig2) const
{
///  Whether or not a given local board has a chance to trig, on each cathode.

  Float_t eff1 = 0.0;
  Float_t eff2 = 0.0;
  GetCellEfficiency(detElemId, localBoard, eff1, eff2);
  trig1 = kTRUE; 
  trig2 = kTRUE;
  if(gRandom->Rndm()>eff1)trig1 = kFALSE;
  if(gRandom->Rndm()>eff2)trig2 = kFALSE;
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFile(const Char_t* filename)
{
///  Reads a file containing the efficiency map.

  TString fileName = gSystem->ExpandPathName(filename);
  if(fileName.EndsWith(".root")){
      ReadHistoBoards(fileName.Data());
      return;
  }

  InitHistos();
  ifstream file(fileName.Data());
  Char_t dat[50];
  if (file.good()){
      file >> dat;
      if(!strcmp(dat,"localBoards"))ReadFileBoards(file);
      else AliWarning("File .dat in wrong format");
      file.close();
  } else {
      AliWarning(Form("Can't read file %s",fileName.Data()));
  }
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFileBoards(ifstream &file)
{
///  Structure of file (.dat) containing local board efficency
  Int_t datInt=0, detEl=0, chamber=0, chCath=0, bin=0;
    Float_t datFloat=0.0;
    Char_t dat[50];

    while (file >> dat) {
	    file >> detEl;
	    chamber = FindChamberIndex(detEl);
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		chCath = fgkNchambers*cath + chamber;
		file >> dat;
		file >> datInt;
		for(Int_t board=1; board<=AliMpConstants::NofLocalBoards(); board++){
		    file >> datFloat;
		    bin = fBoardEfficiency[chCath]->FindBin(board);
		    fBoardEfficiency[chCath]->SetBinContent(bin, datFloat);
		}
	    }
    }
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadHistoBoards(const Char_t *filename)
{
///  Structure of file (.root) containing local board efficency
    TFile *file = new TFile(filename, "read");
    if(!file) {
	AliWarning(Form("Can't read file %s",filename));
	return;
    }
    TString histoName;
    TString cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};
    enum {kAllChEff, kChNonEff, kNumOfHistoTypes};
    TString histoTypeName[2] = {"CountInCh", "NonCountInCh"};

    if(!fCountHistoList) fCountHistoList = new TList();
    else fCountHistoList->Delete();
    if(!fNoCountHistoList) fNoCountHistoList = new TList();
    else fNoCountHistoList->Delete();

    TList *currList[2] = {fCountHistoList, fNoCountHistoList};

    TH1F *histo = 0x0;
    
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      for(Int_t hType=0; hType<kNumOfHistoTypes; hType++){
	histoName = Form("%sChamber%s", cathCode[cath].Data(), histoTypeName[hType].Data());
	histo = (TH1F*)file->Get(histoName.Data());
	currList[hType]->Add(histo);
      }
    }

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t hType=0; hType<kNumOfHistoTypes; hType++){
	  histoName = Form("%sSlat%s%i", cathCode[cath].Data(), histoTypeName[hType].Data(), 11+ch);
	  histo = (TH1F*)file->Get(histoName.Data());
	  currList[hType]->Add(histo);
	}
      }
    }

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t hType=0; hType<kNumOfHistoTypes; hType++){
	  histoName = Form("%sBoard%s%i", cathCode[cath].Data(), histoTypeName[hType].Data(), 11+ch);
	  histo = (TH1F*)file->Get(histoName.Data());
	  currList[hType]->Add(histo);
	}
      }
    }

    InitHistos();
    FillHistosFromList();
}


//_____________________________________________________________________________
void AliMUONTriggerEfficiencyCells::CheckConstants() const
{
/// Check consistence of redefined constants 

  assert(fgkNcathodes == AliMpConstants::NofCathodes());    
  assert(fgkNchambers == AliMpConstants::NofTriggerChambers());    
  assert(fgkNplanes == AliMpConstants::NofTriggerChambers() * fgkNcathodes);    
}


//__________________________________________________________________________
Int_t AliMUONTriggerEfficiencyCells::FindChamberIndex(Int_t detElemId) const
{
///  From detElemId to chamber number

  // Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
  Int_t iChamber = detElemId/100 - 1;
  Int_t chamber = iChamber-AliMpConstants::NofTrackingChambers();
  return chamber;
}


//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::Reset()
{
///  Sets our internal array contents to zero.

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    fBoardEfficiency[chCath] = 0x0;
    fSlatEfficiency[chCath] = 0x0;
  }
}


//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::InitHistos()
{
///  Sets our internal array contents to zero.

  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();
  const Int_t kNslats = 18;
  Int_t chCath=0;
  TString histoName;

  TString cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

  for(Int_t ch=0; ch<fgkNchambers; ch++){
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      chCath = fgkNchambers*cath + ch;
      histoName = Form("%sBoardEffChamber%i", cathCode[cath].Data(), 11+ch);
      fBoardEfficiency[chCath] = new TH1F(histoName.Data(), histoName.Data(), kNumOfBoards, 1-0.5, kNumOfBoards+1.-0.5);
      histoName = Form("%sSlatEffChamber%i", cathCode[cath].Data(), 11+ch);
      fSlatEfficiency[chCath] = new TH1F(histoName.Data(), histoName.Data(), kNslats, 0-0.5, kNslats-0.5);
    }
  }
}


//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::FillHistosFromList()
{
///  Fills internal histos from list.

  Int_t nHistoBins=0;
  TH1F *histoNum = 0x0, *histoDen=0x0, *currHisto = 0x0;
  TString slatName = "Slat", boardName = "Board", histoName;
  Int_t iHistoBoard = -1, iHistoSlat = -1;
  Float_t efficiency, efficiencyError;

  Int_t nentries = fCountHistoList->GetEntries();

  for(Int_t iEntry=0; iEntry<nentries; iEntry++){
    histoNum = (TH1F*)fCountHistoList->At(iEntry);
    histoDen = (TH1F*)fNoCountHistoList->At(iEntry);

    if(!histoNum) {
      AliWarning("Histogram not found in fCountHistoList. Skip to next");
      continue;
    }
    if(!histoDen) {
      AliWarning("Histogram not found in fNoCountHistoList. Skip to next");
      continue;
    }

    histoName = histoNum->GetName();
    nHistoBins = histoNum->GetNbinsX();

    if(histoName.Contains(boardName)){
      iHistoBoard++;
      currHisto = fBoardEfficiency[iHistoBoard];
    }
    else if(histoName.Contains(slatName)){
      iHistoSlat++;
      currHisto = fSlatEfficiency[iHistoSlat];
    }
    else continue;

    for(Int_t iBin=1; iBin<=nHistoBins; iBin++){
      CalculateEfficiency((Int_t)histoNum->GetBinContent(iBin), (Int_t)histoNum->GetBinContent(iBin) + (Int_t)histoDen->GetBinContent(iBin), efficiency, efficiencyError, kFALSE);

      currHisto->SetBinContent(iBin, efficiency);
      currHisto->SetBinError(iBin, efficiencyError);
    }
  }
}


//_____________________________________________________________________________
void AliMUONTriggerEfficiencyCells::CalculateEfficiency(Int_t trigger44, Int_t trigger34,
							Float_t &efficiency, Float_t &error,
							Bool_t failuresAsInput)
{
    //
    /// Returns the efficiency.
    //

    efficiency=-9.;
    error=0.;
    if(trigger34>0){
	efficiency=(Double_t)trigger44/((Double_t)trigger34);
	if(failuresAsInput)efficiency=1.-(Double_t)trigger44/((Double_t)trigger34);
    }
    Double_t q = TMath::Abs(1-efficiency);
    if(efficiency<0)error=0.0;
    else error = TMath::Sqrt(efficiency*q/((Double_t)trigger34));
}


//_____________________________________________________________________________
void AliMUONTriggerEfficiencyCells::CheckFiredStrips(const Char_t* cdbStorage,
						     Int_t runNumber)
{
  //
  /// Check for fired strips participating to efficiency
  /// calculation (when available).
  /// Strips inside a local board should be quite homogeneously hit
  /// If not, this could be a problem of electronics (i.e. ADULT board off).
  //

  if(!fFiredStrips) {
    AliWarning("List of fired pads not present. Check not performable.");
    return;
  }

  GetListsForCheck(cdbStorage, runNumber);

  TString histoName;

  // Check fired pads (when available)
  if(fFiredFitHistoList){
    TH1F *histo1D = 0x0;
    TF1 *fitFunc = 0x0;
    TCanvas *histoFiredCan[20];
    Int_t nEntries = fFiredFitHistoList->GetEntries();
    for(Int_t iEntry=0; iEntry<nEntries; iEntry++){
      histo1D = (TH1F*)fFiredFitHistoList->At(iEntry);
      printf("Problems found in %s\n", histo1D->GetTitle());
    }
    Int_t nPrintCan = nEntries;
    if(nPrintCan>20) {
      AliWarning("Too many boards with problems: only 20 will be shown");
      nPrintCan = 20;
    }
    for(Int_t iCan=0; iCan<nPrintCan; iCan++){
      histo1D = (TH1F*)fFiredFitHistoList->At(iCan);
      histoName = histo1D->GetName();
      histoName.Append("Can");
      histoFiredCan[iCan] = new TCanvas(histoName.Data(), histoName.Data(), 100+10*iCan, 10*iCan, 700, 700);
      histoFiredCan[iCan]->SetRightMargin(0.14);
      histoFiredCan[iCan]->SetLeftMargin(0.12);
      histo1D->Draw("E");
      fitFunc = histo1D->GetFunction("pol0");
      fitFunc->SetLineColor(2);
      fitFunc->Draw("same");
    }
    if(nEntries==0){
      printf("\nAll local boards seem ok!!\n\n");
    }
  }
}


//_____________________________________________________________________________
void AliMUONTriggerEfficiencyCells::DisplayEfficiency(Bool_t perSlat,
						      const Char_t* cdbStorage,
						      Int_t runNumber)
{
  //
  /// Display calculated efficiency.
  //

  Bool_t isInitSlat = kFALSE, isInitBoard = kFALSE;
  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    if(fBoardEfficiency[chCath]->GetEntries()>0) isInitBoard = kTRUE;
    if(fSlatEfficiency[chCath]->GetEntries()>0) isInitSlat = kTRUE;
  }

  if(!isInitBoard){
    printf("Trigger efficiency not initialized per board.\nDisplay not yet implemented.\n");
    return;
  }
  if(!isInitSlat && perSlat){
    printf("Trigger efficiency not initialized for slat.\nPlease try option kFALSE.\n");
    return;
  }
  
  GetListsForCheck(cdbStorage, runNumber);

  //const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

  TH2F *histo = 0x0;
  TString histoName, histoTitle;

  // Plot fired strips (when available)
  if(fFiredDisplayHistoList){
    TCanvas *displayFiredCan[fgkNplanes];
    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      histo = (TH2F*)fFiredDisplayHistoList->At(chCath);
      histoName = Form("%sCan", histo->GetName());
      histoTitle = Form("%s", histo->GetTitle());
      displayFiredCan[chCath] = new TCanvas(histoName.Data(), histoTitle.Data(), 100+10*chCath, 10*chCath, 700, 700);
      displayFiredCan[chCath]->SetRightMargin(0.14);
      displayFiredCan[chCath]->SetLeftMargin(0.12);
      histo->GetYaxis()->SetTitleOffset(1.4);
      histo->SetStats(kFALSE);
      histo->Draw("COLZ");
    }
  }

  // Plot efficiency
  if(fDisplayHistoList){
    TCanvas *can[fgkNplanes];
    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      Int_t currChCath = chCath;
      if(perSlat==kTRUE) currChCath += fgkNplanes;
      histo = (TH2F*)fDisplayHistoList->At(currChCath);
      histoName = Form("%sCan", histo->GetName());
      histoTitle = Form("%s", histo->GetTitle());
      can[chCath] = new TCanvas(histoName.Data(), histoTitle.Data(), 100+10*chCath, 10*chCath, 700, 700);
      can[chCath]->SetRightMargin(0.14);
      can[chCath]->SetLeftMargin(0.12);
      histo->GetZaxis()->SetRangeUser(0.,1.);
      histo->GetYaxis()->SetTitleOffset(1.4);
      histo->SetStats(kFALSE);
      histo->Draw("COLZ");
      if(perSlat==kTRUE) continue;
      histo = (TH2F*)fBoardLabelList->At(currChCath);
      histo->Draw("textsame");
    }
  }
}


//__________________________________________________________________________
Bool_t AliMUONTriggerEfficiencyCells::GetListsForCheck(const Char_t* cdbStorage,
						       Int_t runNumber)
{
  //
  /// Getting histograms for efficiency, 
  /// map of fired strips entering efficiency calculations,
  /// fits for checking switched-off elements in chambers.
  //
  const Float_t kChi2RedMax = 1.5;
  const Float_t kDummyFired = 1e-5;

  if(fDisplayHistoList || fBoardLabelList || fFiredFitHistoList || fFiredDisplayHistoList) return kTRUE;

  if(!fDisplayHistoList) fDisplayHistoList = new TList();
  if(!fBoardLabelList) fBoardLabelList = new TList(); 
  if(!fFiredFitHistoList && fFiredStrips) fFiredFitHistoList = new TList();
  if(!fFiredDisplayHistoList && fFiredStrips) fFiredDisplayHistoList = new TList();

  AliCDBManager::Instance()->SetDefaultStorage(cdbStorage);
  AliCDBManager::Instance()->SetRun(runNumber);

  TH3F* padFired = 0x0;
  TH2F* displayHisto = 0x0;
  TH1F* histoFired[fgkNplanes][234];
  TF1* fitFunc = 0x0;
  Bool_t isStripOffInBoard[fgkNplanes][234];
  Int_t bin=0;
  TString histoName, histoTitle;
  TString cathName[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

  AliMUONTriggerDisplay triggerDisplay;
  // Book histos
  for(Int_t iCath=0; iCath<fgkNcathodes; iCath++){
    if(fFiredStrips) padFired = (TH3F*)fFiredStrips->At(iCath);
    for(Int_t iCh=0; iCh<fgkNchambers; iCh++){
      Int_t chCath = fgkNchambers*iCath + iCh;
      Int_t currCh = 11 + iCh;
      histoName = Form("%sChamber%i", cathName[iCath].Data(), currCh);
      histoTitle = Form("Chamber %i: efficiency %s", currCh, cathName[iCath].Data());
      displayHisto = 
	(TH2F*)triggerDisplay.GetDisplayHistogram(fBoardEfficiency[chCath], histoName,
						  AliMUONTriggerDisplay::kDisplayBoards,
						  iCath,currCh,histoTitle,
						  AliMUONTriggerDisplay::kShowZeroes);
      fDisplayHistoList->Add(displayHisto);

      histoName = Form("labels%sChamber%i", cathName[iCath].Data(), currCh);
      displayHisto = 
	(TH2F*)triggerDisplay.GetBoardNumberHisto(histoName,currCh);
      fBoardLabelList->Add(displayHisto);

      if(!fFiredStrips) continue;
      histoName = Form("firedPads%sChamber%i", cathName[iCath].Data(), currCh);
      histoTitle = Form("Chamber %i: Fired pads %s", currCh, cathName[iCath].Data());
      bin = padFired->GetXaxis()->FindBin(currCh);
      padFired->GetXaxis()->SetRange(bin,bin);
      displayHisto = 
	(TH2F*)triggerDisplay.GetDisplayHistogram(padFired->Project3D("zy"),histoName,
						  AliMUONTriggerDisplay::kDisplayStrips,
						  iCath,currCh,histoTitle);
      fFiredDisplayHistoList->Add(displayHisto);
      
      for(Int_t ib=0; ib<AliMpConstants::NofLocalBoards(); ib++){
	histoName = Form("%sChamber%iBoard%i", cathName[iCath].Data(), currCh, ib+1);
	histoTitle = Form("Chamber %i: fired pads %s board = %i", currCh, cathName[iCath].Data(), ib+1);
	histoFired[chCath][ib] = new TH1F(histoName.Data(), histoTitle.Data(), 16, -0.5, 15.5);
	histoFired[chCath][ib]->SetXTitle("board");
	isStripOffInBoard[chCath][ib] = kFALSE;
      } // loop on board
    } // loop on chamber
  } // loop on cathode

  for(Int_t iCath=0; iCath<fgkNcathodes; iCath++){
    for(Int_t iCh=0; iCh<fgkNchambers; iCh++){
      Int_t chCath = fgkNchambers*iCath + iCh;
      Int_t currCh = 11+iCh;
      histoName = Form("%sChamber%iSlatEff", cathName[iCath].Data(), currCh);
      histoTitle = Form("Chamber %i: efficiency %s per slat", currCh, cathName[iCath].Data());
      displayHisto = 
	(TH2F*)triggerDisplay.GetDisplayHistogram(fSlatEfficiency[chCath], histoName,
						  AliMUONTriggerDisplay::kDisplaySlats,
						  iCath,currCh,histoTitle);
      fDisplayHistoList->Add(displayHisto);
    }
  }

  if(!fFiredStrips) return kTRUE;

  // Check fired pads (when available)
  for(Int_t iLoc = 0; iLoc < AliMpConstants::NofLocalBoards(); iLoc++) {  
    Int_t iBoard = iLoc+1;
    for(Int_t iCh=0; iCh<AliMpConstants::NofChambers(); iCh++){
      Int_t iChamber = iCh + AliMpConstants::NofTrackingChambers();

      Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(iBoard, iChamber);

      if (!detElemId) continue;

      AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(iBoard, kFALSE);

      // skip copy cards
      if( !localBoard->IsNotified()) 
	continue;
      
      for(Int_t iCath=0; iCath<AliMpConstants::NofCathodes(); iCath++){
	Int_t chCath = fgkNchambers*iCath + iCh;
	// loop over strips
	const AliMpVSegmentation* seg = 
	  AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(iCath));
	Int_t nStrips=0;
	for (Int_t iStrip = 0; iStrip < 16; ++iStrip) {
	  AliMpPad pad = seg->PadByLocation(iBoard,iStrip,kFALSE);
	  if (!pad.IsValid()) continue;
	  nStrips++;
	  padFired = (TH3F*)fFiredStrips->At(iCath);
	  Float_t nFired = padFired->GetBinContent(11+iCh, iBoard, iStrip);
	  if(nFired==0.) nFired = kDummyFired;
	  histoFired[chCath][iLoc]->Fill(iStrip, nFired);
	}

	histoFired[chCath][iLoc]->Fit("pol0","Q0R","",0., (Float_t)nStrips-1.);
	fitFunc = histoFired[chCath][iLoc]->GetFunction("pol0");
	Float_t chi2 = fitFunc->GetChisquare();
	Float_t ndf = (Float_t)fitFunc->GetNDF();
	Float_t reducedChi2 = chi2/ndf;
	if(reducedChi2>kChi2RedMax) {
	  isStripOffInBoard[chCath][iLoc] = kTRUE;
	  fFiredFitHistoList->Add(histoFired[chCath][iLoc]);
	}
      } // loop on cathodes
    } // loop on chambers
  } // loop on local boards 

  return kTRUE;
}

//__________________________________________________________________________
Bool_t AliMUONTriggerEfficiencyCells::SumRunEfficiency(const AliMUONTriggerEfficiencyCells &other)
{
///  Sums results from different runs and gives the efficiency
  if(!fCountHistoList || !fNoCountHistoList) {
    AliWarning("Histograms for efficiency calculations not implemented in object");
    return kFALSE;
  }
  if(!other.fCountHistoList || !other.fNoCountHistoList) {
    AliWarning("Histograms for efficiency calculations not implemented in object passed as argument");
    return kFALSE;
  }

  Int_t nentries = fCountHistoList->GetEntries();
  TH1F *currNum = 0x0, *currDen = 0x0, *otherNum = 0x0, *otherDen = 0x0;

  for(Int_t iEntry=0; iEntry<nentries; iEntry++){
    currNum = (TH1F*)fCountHistoList->At(iEntry);
    currDen = (TH1F*)fNoCountHistoList->At(iEntry);
    otherNum = (TH1F*)other.fCountHistoList->At(iEntry);
    otherDen = (TH1F*)other.fNoCountHistoList->At(iEntry);
    currNum->Add(otherNum);
    currDen->Add(otherDen);
  }

  FillHistosFromList();

  if(!fFiredStrips) {
    AliWarning("Histograms for fired region check not implemented in object");
    return kFALSE;
  }
  if(!other.fFiredStrips) {
    AliWarning("Histograms for fired region check not implemented in object passed as argument");
    return kFALSE;
  }
  
  TH3F *currFired = 0x0, *otherFired = 0x0;
  nentries = fFiredStrips->GetEntries();
  for(Int_t iEntry=0; iEntry<nentries; iEntry++){
    currFired  = (TH3F*)fFiredStrips->At(iEntry);
    otherFired = (TH3F*)fFiredStrips->At(iEntry);
    currFired->Add(otherFired);
  }
    
  return kTRUE;
}
