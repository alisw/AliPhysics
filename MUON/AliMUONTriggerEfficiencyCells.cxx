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

#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMpConstants.h"
#include "AliMpDEManager.h"

// Classes for display
#include "AliMUONGeometryTransformer.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDDL.h"
#include "AliMpTriggerCrate.h"
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
/// Efficiency is stored per cathode on local boards, or, alternatively,
/// on "cells" of a given size.
///
/// The main method of this class is IsTriggered().
///
/// $ALICE_ROOT/MUON/data/TriggerChamberefficiencyCells.dat contains efficiency 
/// for each chamber (i.e. DetElement). 
///
/// In the case of local boards, efficiency is stored from left to right
/// per increasing board number (from 1 to 234)
///
/// Otherwise, he efficiency cells goes from right to left and 
/// from bottom to top of the chamber, namely, the efficiencies tabulated in the 
/// file refers to the following reference frame:
///
/// <pre>
/// x
/// <----------------------------------|
///                                    |
///    ---------------------------     |
///   | 0.97 | 0.97 | 0.97 | 0.97 |    |
///    ---------------------------     |
///   | 0.97 | 0.97 | 0.97 | 0.97 |    |
///    ---------------------------     |
///   | 0.97 | 0.97 | 0.97 | 0.97 |    |
///    ---------------------------     |
///                                    |
///                                   \/ y
/// </pre>
///
///  In both cases, the file can be edited in order to change efficiency
///  in a chosen local board/region of the chamber.
///
///
/// But please note that this object is also available from the CDB 
/// (generated using the MUONCDB.C macro)
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
    for(Int_t slat=0; slat<fgkNslats; slat++){
      fCellContent[chCath][slat] = other.fCellContent[chCath][slat];
    }
    fCellSize[chCath] = other.fCellSize[chCath];
    fCellNumber[chCath] = other.fCellNumber[chCath];
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
    for(Int_t slat=0; slat<fgkNslats; slat++){
      fCellContent[chCath][slat] = other.fCellContent[chCath][slat];
    }
    fCellSize[chCath] = other.fCellSize[chCath];
    fCellNumber[chCath] = other.fCellNumber[chCath];
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
void AliMUONTriggerEfficiencyCells::GetCellEfficiency(Int_t detElemId, Float_t x, Float_t y, Float_t &eff1, Float_t &eff2) const
{
///  Get the efficiencies of the 2 cathodes at a given location (x,y)

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t slat = FindSlatIndex(detElemId);
  TArrayI cell = CellByCoord(detElemId,x,y);
  eff1 = 0.0;
  eff2 = 0.0;
  if(cell.At(0)>=0 && cell.At(1)>=0)
  {
    eff1 = fCellContent[chamber][slat][cell.At(0)][cell.At(1)];
    eff2 = fCellContent[fgkNchambers+chamber][slat][cell.At(0)][cell.At(1)];
  }
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
AliMUONTriggerEfficiencyCells::IsTriggered(Int_t detElemId, Float_t x, Float_t y, Bool_t &trig1, Bool_t &trig2) const
{
///  Whether or not a given location (x,y) has a chance to trig, on each cathode.

  Float_t eff1 = 0.0;
  Float_t eff2 = 0.0;
  GetCellEfficiency(detElemId, x, y, eff1, eff2);
  trig1 = kTRUE; 
  trig2 = kTRUE;
  if(gRandom->Rndm()>eff1)trig1 = kFALSE;
  if(gRandom->Rndm()>eff2)trig2 = kFALSE;
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
TArrayI AliMUONTriggerEfficiencyCells::CellByCoord(Int_t detElemId, Float_t x, Float_t y) const
{
///  Get the efficiencies at a given location.

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t slat = FindSlatIndex(detElemId);
  Int_t cell[fgkNcathodes]={-1,-1};
  Float_t maxX = fCellSize[chamber][slat]*((Float_t)fCellNumber[chamber][slat]);
  Float_t maxY = fCellSize[fgkNchambers+chamber][slat]*((Float_t)fCellNumber[fgkNchambers+chamber][slat]);
  if(x>=0 & x<maxX & y>=0 & y<maxY)
  {
    cell[0] = (Int_t)(x/fCellSize[chamber][slat]);
    cell[1] = (Int_t)(y/fCellSize[fgkNchambers+chamber][slat]);
  }
  return TArrayI(fgkNcathodes,cell);
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
      else ReadFileXY(file);
      file.close();
  } else {
      AliWarning(Form("Can't read file %s",fileName.Data()));
  }
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFileXY(ifstream &file)
{
///  Structure of file (.dat) containing geometrical efficency
    Int_t datInt=0, detEl=0, chamber=0, rpc=0, chCath=0;
    Float_t datFloat=0.0;
    Char_t dat[50];

    while (file >> dat) {
	file >> detEl;
	chamber = FindChamberIndex(detEl);
	rpc = FindSlatIndex(detEl);
	file >> dat;
	for(Int_t i=0; i<fgkNcathodes; i++){
	    chCath = fgkNchambers*i + chamber;
	    file >> datInt;
	    fCellNumber[chCath][rpc] = datInt;
	    file >> dat;
	}
	for(Int_t i=0; i<fgkNcathodes; i++){
	    chCath = fgkNchambers*i + chamber;
	    file >> datFloat;
	    fCellSize[chCath][rpc] = datFloat;
	    if(i==0)file >> dat;
	}
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    chCath = fgkNchambers*cath + chamber;
	    file >> dat;
	    file >> datInt;
	    for(Int_t iy=0; iy<fCellNumber[fgkNchambers+chamber][rpc]; iy++){
		for(Int_t ix=0; ix<fCellNumber[chamber][rpc]; ix++){
		    file >> datFloat;
		    fCellContent[chCath][rpc][ix][iy] = datFloat;
		}
	    }
	}
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
	    //rpc = FindSlatIndex(detEl);
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
    Char_t histoName[30];
    Char_t *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

    for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    sprintf(histoName, "%sBoardEffChamber%i", cathCode[cath], 11+ch);
	    if(!(TH1F *)file->Get(histoName)) {
		AliWarning(Form("Can't find histo %s in file %s",histoName, filename));
		continue;
	    }
	    Int_t chCath = fgkNchambers*cath + ch;
	    fBoardEfficiency[chCath] = (TH1F *)file->Get(histoName);
	}
    }
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
  Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
  Int_t chamber = iChamber-AliMpConstants::NofTrackingChambers();
  return chamber;
}


//__________________________________________________________________________
Int_t AliMUONTriggerEfficiencyCells::FindSlatIndex(Int_t detElemId) const
{
///  From detElemId to slat index.
  Int_t slat = detElemId%100;
  return slat;
}


//__________________________________________________________________________
TVector2 AliMUONTriggerEfficiencyCells::ChangeReferenceFrame(Float_t x, Float_t y, Float_t x0, Float_t y0)
{
/// (x0,y0) position of the local reference frame (center of the chamber)

    Float_t x1 = x0-x;//reflection of axis
    Float_t y1 = y+y0;
    return TVector2(x1,y1);
}

//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::Reset()
{
///  Sets our internal array contents to zero.

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fCellSize[chCath].Set(fgkNslats);
	fCellNumber[chCath].Set(fgkNslats);
	for(Int_t slat=0; slat<fgkNslats; slat++){
	    fCellContent[chCath][slat].ResizeTo(fgkNcells,fgkNcells);
	}
    }

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fCellSize[chCath].Reset();
	fCellNumber[chCath].Reset();
	for(Int_t slat=0; slat<fgkNslats; slat++){
	    fCellContent[chCath][slat].Zero();
	}
    }

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      fBoardEfficiency[chCath] = 0x0;
      fSlatEfficiency[chCath] = 0x0;
    }

    if(!AliMpCDB::LoadDDLStore()) // Load DDL store from OCDB
      AliFatal("Could not access mapping OCDB");
}


//__________________________________________________________________________
void
AliMUONTriggerEfficiencyCells::InitHistos()
{
///  Sets our internal array contents to zero.

  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();
  Int_t chCath=0;
  Char_t histoName[40];

  Char_t *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

  for(Int_t ch=0; ch<fgkNchambers; ch++){
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      chCath = fgkNchambers*cath + ch;
      sprintf(histoName, "%sBoardEffChamber%i", cathCode[cath], 11+ch);
      fBoardEfficiency[chCath] = new TH1F(histoName, histoName, kNumOfBoards, 1-0.5, kNumOfBoards+1.-0.5);
      sprintf(histoName, "%sSlatEffChamber%i", cathCode[cath], 11+ch);
      fSlatEfficiency[chCath] = new TH1F(histoName, histoName, fgkNslats, 0-0.5, fgkNslats-0.5);
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
void AliMUONTriggerEfficiencyCells::CheckFiredStrips(const Char_t* geoFilename)
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

  GetListsForCheck(geoFilename);

  Char_t histoName[40], histoTitle[90];

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
      histoFiredCan[iCan] = new TCanvas(histoName, histoTitle, 100+10*iCan, 10*iCan, 700, 700);
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
void AliMUONTriggerEfficiencyCells::DisplayEfficiency(Bool_t perSlat, const Char_t* geoFilename)
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
  
  GetListsForCheck(geoFilename);

  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

  TH2F *histo = 0x0;
  Char_t histoName[40], histoTitle[90];
  TPaveLabel *boardLabel = 0x0;

  // Plot fired strips (when available)
  if(fFiredDisplayHistoList){
    TCanvas *displayFiredCan[fgkNplanes];
    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      histo = (TH2F*)fFiredDisplayHistoList->At(chCath);
      sprintf(histoName, "%sCan", histo->GetName());
      sprintf(histoTitle, "%s", histo->GetTitle());
      displayFiredCan[chCath] = new TCanvas(histoName, histoTitle, 100+10*chCath, 10*chCath, 700, 700);
      displayFiredCan[chCath]->SetRightMargin(0.14);
      displayFiredCan[chCath]->SetLeftMargin(0.12);
      histo->GetYaxis()->SetTitleOffset(1.4);
      histo->SetStats(kFALSE);
      histo->Draw("COLZ");
      for (Int_t board = 0; board < kNumOfBoards; board++) {
	Int_t currLabel = chCath * kNumOfBoards + board;
	boardLabel = (TPaveLabel*)fBoardLabelList->At(currLabel);
	boardLabel->Draw("same");
      }
    }
  }

  // Plot efficiency
  if(fDisplayHistoList){
    TCanvas *can[fgkNplanes];
    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      Int_t currChCath = chCath;
      if(perSlat==kTRUE) currChCath += fgkNplanes;
      histo = (TH2F*)fDisplayHistoList->At(currChCath);
      sprintf(histoName, "%sCan", histo->GetName());
      sprintf(histoTitle, "%s", histo->GetTitle());
      can[chCath] = new TCanvas(histoName, histoTitle, 100+10*chCath, 10*chCath, 700, 700);
      can[chCath]->SetRightMargin(0.14);
      can[chCath]->SetLeftMargin(0.12);
      histo->GetZaxis()->SetRangeUser(0.,1.);
      histo->GetYaxis()->SetTitleOffset(1.4);
      histo->SetStats(kFALSE);
      histo->Draw("COLZ");
      for (Int_t board = 0; board < kNumOfBoards; board++) {
	Int_t currLabel = chCath * kNumOfBoards + board;
	if(perSlat==kTRUE) currLabel += kNumOfBoards * fgkNplanes;
	boardLabel = (TPaveLabel*)fBoardLabelList->At(currLabel);
	boardLabel->Draw("same");
      }
    }
  }
}


//__________________________________________________________________________
Bool_t AliMUONTriggerEfficiencyCells::GetListsForCheck(const Char_t* geoFilename)
{
  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();
  const Float_t kChi2RedMax = 1.5;
  const Float_t kDummyFired = 1e-5;

  if(fDisplayHistoList || fBoardLabelList || fFiredFitHistoList || fFiredDisplayHistoList) return kTRUE;

  if(!fDisplayHistoList) fDisplayHistoList = new TList();
  if(!fBoardLabelList) fBoardLabelList = new TList(); 
  if(!fFiredFitHistoList && fFiredStrips) fFiredFitHistoList = new TList();
  if(!fFiredDisplayHistoList && fFiredStrips) fFiredDisplayHistoList = new TList();

  AliMUONGeometryTransformer *transform = new AliMUONGeometryTransformer();
  transform->LoadGeometryData(geoFilename);

  if(!AliMpCDB::LoadDDLStore()) // Load DDL store from OCDB
    AliFatal("Could not access mapping OCDB");

  AliMpDDLStore *ddlStore = AliMpDDLStore::Instance();
  Int_t line, slat;
  Float_t xLocal1=0., yLocal1=0., xLocal2=0., yLocal2=0.;
  Float_t xg1, yg1, zg1, xg2, yg2, zg2;
  Float_t xWidth=0., yWidth=0.;
  Float_t x1Label=0., x2Label=0., y1Label=0., y2Label=0.;
  Float_t x1LabelSlat=0., x2LabelSlat=0., y1LabelSlat=0., y2LabelSlat=0.;
  Int_t x1=0, y1=0, x2=0, y2=0, localId;

  gStyle->SetPalette(1);

  Char_t *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

  Float_t boardsX = 280.00;  // cm
  Float_t boardsY = 335.00;  // cm

  // Check fired pads (when available)  
  Int_t maxY[fgkNcathodes] = {64,1};
  TH3F *padFired[fgkNplanes];
  TH1F *histoFired[fgkNplanes][234];
  TH2F *histoFiredDisplay[fgkNplanes];
  TF1 *fitFunc = 0x0;
  Bool_t isStripOffInBoard[fgkNplanes][234];

  TH2F *histo[fgkNplanes];
  TH2F *histoSlat[fgkNplanes];
  TPaveLabel *boardLabel[fgkNplanes][234];
  TPaveLabel *boardLabelSlat[fgkNplanes][234];
  assert(kNumOfBoards==234);

  Char_t histoName[40], histoTitle[90], labelTxt[5], labelSlatTxt[5];

  Float_t efficiency, efficiencyError, efficiencySlat, efficiencySlatError;

  // Book histos
  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    padFired[cath] = 0x0;
    if(fFiredStrips) padFired[cath] = (TH3F*)fFiredStrips->At(cath);
    for(Int_t ch=0; ch<fgkNchambers; ch++){
      Int_t chCath = fgkNchambers*cath + ch;
      sprintf(histoName, "%sChamber%i", cathCode[cath], 11+ch);
      sprintf(histoTitle, "Chamber %i: efficiency %s", 11+ch, cathCode[cath]);
      histo[chCath] = new TH2F(histoName, histoTitle, (Int_t)boardsX, -boardsX, boardsX, (Int_t)boardsY, -boardsY, boardsY);
      histo[chCath]->SetXTitle("X (cm)");
      histo[chCath]->SetYTitle("Y (cm)");
      fDisplayHistoList->Add(histo[chCath]);

      if(!fFiredStrips) continue;
      sprintf(histoName, "firedPads%sChamber%i", cathCode[cath], 11+ch);
      sprintf(histoTitle, "Chamber %i: Fired pads %s", 11+ch, cathCode[cath]);
      histoFiredDisplay[chCath] = new TH2F(histoName, histoTitle, (Int_t)boardsX, -boardsX, boardsX, (Int_t)boardsY, -boardsY, boardsY);
      histoFiredDisplay[chCath]->SetXTitle("X (cm)");
      histoFiredDisplay[chCath]->SetYTitle("Y (cm)");
      fFiredDisplayHistoList->Add(histoFiredDisplay[chCath]);
      
      for(Int_t ib=0; ib<kNumOfBoards; ib++){
	sprintf(histoName, "%sChamber%iBoard%i", cathCode[cath], 11+ch, ib+1);
	sprintf(histoTitle, "Chamber %i: fired pads %s board = %i", 11+ch, cathCode[cath], ib+1);
	histoFired[chCath][ib] = new TH1F(histoName, histoTitle, 16, -0.5, 15.5);
	histoFired[chCath][ib]->SetXTitle("board");
	isStripOffInBoard[chCath][ib] = kFALSE;
      } // loop on board
    } // loop on chamber
  } // loop on cathode

  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    for(Int_t ch=0; ch<fgkNchambers; ch++){
      Int_t chCath = fgkNchambers*cath + ch;
      sprintf(histoName, "%sChamber%iSlatEff", cathCode[cath], 11+ch);
      sprintf(histoTitle, "Chamber %i: efficiency %s per slat", 11+ch, cathCode[cath]);
      histoSlat[chCath] = new TH2F(histoName, histoTitle, (Int_t)boardsX, -boardsX, boardsX, (Int_t)boardsY, -boardsY, boardsY);
      histoSlat[chCath]->SetXTitle("X (cm)");
      histoSlat[chCath]->SetYTitle("Y (cm)");
      fDisplayHistoList->Add(histoSlat[chCath]);
    }
  }

  // loop over the trigger DDL (Right: 20, Left: 21)
  for (Int_t iDDL = 20; iDDL <= 21; ++iDDL) {
    AliMpDDL* ddl = ddlStore->GetDDL(iDDL);
    Int_t nCrate = ddl->GetNofTriggerCrates();
    // loop over the number of crates in DDL
    for (Int_t index = 0; index < nCrate; ++index) {
      // get crate object
      AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, index);
      Int_t nLocal = crate->GetNofLocalBoards();
      for (Int_t iLocal = 0; iLocal < nLocal; ++iLocal) {
	// get local board Id from crate object
	localId = crate->GetLocalBoardId(iLocal);
	// get local board object
	AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(localId);

	if (!localBoard->IsNotified()) continue;

	// get detection element connected to this board
	for (Int_t ch = 0; ch < fgkNchambers; ++ch) {
	  Int_t iCh = ch + AliMpConstants::NofTrackingChambers();
	  Int_t detElemId = ddlStore->GetDEfromLocalBoard(localId, iCh);

	  if (!detElemId) continue;

	  // get segmentation
	  for (Int_t cath = 0; cath < fgkNcathodes; ++cath) {
	    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()
	      ->GetMpSegmentation(detElemId,
				  AliMp::GetCathodType(cath));

	    Int_t chCath = fgkNchambers*cath + ch;
	    slat = detElemId%100;
	    Int_t nStrips=0;

	    // loop over strips
	    for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) {
	      // get pad from electronics
	      AliMpPad pad = seg->PadByLocation(AliMpIntPair(localId,
							       ibitxy),kFALSE);
		
	      if (!pad.IsValid()) continue;
	      if(cath==0){ // Geometry info from bending plane only
		if(ibitxy==0) {
		  xLocal1 = pad.Position().X();
		  yLocal1 = pad.Position().Y();
		  xWidth = pad.Dimensions().X();
		  yWidth = pad.Dimensions().Y();
		}
		xLocal2 = pad.Position().X();
		yLocal2 = pad.Position().Y();
	      }
	      
	      // Check fired pads (when available)
	      if(padFired[cath]) {
		Int_t padX = pad.GetIndices().GetFirst();
		Int_t padY = pad.GetIndices().GetSecond();
		Int_t currPair = padX*maxY[cath] + padY;
		nStrips++;
		//printf("cath %i board = %i  (%2i, %2i) -> %i\n",
		//cath, localId, padX, padY, currPair);
		Int_t chBin = padFired[cath]->GetXaxis()->FindBin(ch);
		Int_t slatBin = padFired[cath]->GetYaxis()->FindBin(slat);
		Int_t pairBin = padFired[cath]->GetZaxis()->FindBin(currPair);
		Float_t nFired = padFired[cath]->GetBinContent(chBin, slatBin, pairBin);;

		Float_t dimX = pad.Dimensions().X();
		Float_t dimY = pad.Dimensions().Y();

		Float_t stripX1 = pad.Position().X();
		Float_t stripY1 = pad.Position().Y();
		Float_t stripX2 = pad.Position().X();
		Float_t stripY2 = pad.Position().Y();

		transform->Local2Global(detElemId, stripX1, stripY1, 0, xg1, yg1, zg1);
		transform->Local2Global(detElemId, stripX2, stripY2, 0, xg2, yg2, zg2);

		Float_t x1Float = TMath::Min(xg1,xg2) - dimX;
		Float_t y1Float = TMath::Min(yg1,yg2) - dimY;
		Float_t x2Float = TMath::Max(xg1,xg2) + dimX;
		Float_t y2Float = TMath::Max(yg1,yg2) + dimY;  

		Int_t x1Int = histoFiredDisplay[chCath]->GetXaxis()->FindBin(x1Float)+1;
		Int_t y1Int = histoFiredDisplay[chCath]->GetYaxis()->FindBin(y1Float)+1;
		Int_t x2Int = histoFiredDisplay[chCath]->GetXaxis()->FindBin(x2Float)-1;
		Int_t y2Int = histoFiredDisplay[chCath]->GetYaxis()->FindBin(y2Float)-1;

		for(Int_t binX=x1Int; binX<=x2Int; binX++){
		  for(Int_t binY=y1Int; binY<=y2Int; binY++){
		    histoFiredDisplay[chCath]->SetBinContent(binX, binY, nFired);
		  }
		}

		if(nFired==0.) nFired = kDummyFired;
		histoFired[chCath][localId-1]->Fill(ibitxy, nFired);
	      }
	    } // loop on strips

	    if(cath==0){ // Geometry info from bending plane only
	      transform->Local2Global(detElemId, xLocal1, yLocal1, 0, xg1, yg1, zg1);
	      transform->Local2Global(detElemId, xLocal2, yLocal2, 0, xg2, yg2, zg2);

	      // Per board
	      x1Label = TMath::Min(xg1,xg2) - xWidth;
	      y1Label = TMath::Min(yg1,yg2) - yWidth;
	      x2Label = TMath::Max(xg1,xg2) + xWidth;
	      y2Label = TMath::Max(yg1,yg2) + yWidth;

	      x1 = histo[ch]->GetXaxis()->FindBin(x1Label)+1;
	      y1 = histo[ch]->GetYaxis()->FindBin(y1Label)+1;
	      x2 = histo[ch]->GetXaxis()->FindBin(x2Label)-1;
	      y2 = histo[ch]->GetYaxis()->FindBin(y2Label)-1;

	      sprintf(labelTxt,"%3d", localId);

	      // Per slat
	      line = localBoard->GetPosition().GetFirst();
	      x1LabelSlat = 140.;
	      x2LabelSlat = x1LabelSlat + 40.;
	      y1LabelSlat = -285. + ((Float_t)(line - 1)) * 68;
	      y2LabelSlat = y1LabelSlat + 34.;
	      if(localId>kNumOfBoards/2){
		x1LabelSlat = -x2LabelSlat;
		x2LabelSlat = x1LabelSlat + 40.;
	      }
	      sprintf(labelSlatTxt,"%2d", slat);
	    }

	    boardLabel[chCath][localId-1] = new TPaveLabel(x1Label, y1Label, x2Label, y2Label, labelTxt);
	    boardLabel[chCath][localId-1]->SetFillStyle(0);
	    boardLabel[chCath][localId-1]->SetBorderSize(0);

	    boardLabelSlat[chCath][localId-1] = new TPaveLabel(x1LabelSlat, y1LabelSlat, x2LabelSlat, y2LabelSlat, labelSlatTxt);
	    boardLabelSlat[chCath][localId-1]->SetFillStyle(0);
	    boardLabelSlat[chCath][localId-1]->SetBorderSize(0);

	    Int_t histoBin = localId;
	    efficiency = fBoardEfficiency[chCath]->GetBinContent(histoBin);
	    efficiencyError = fBoardEfficiency[chCath]->GetBinError(histoBin);

	    histoBin = slat+1;
	    efficiencySlat = fSlatEfficiency[chCath]->GetBinContent(histoBin);
	    efficiencySlatError = fSlatEfficiency[chCath]->GetBinError(histoBin);

	    for(Int_t binX=x1; binX<=x2; binX++){
	      for(Int_t binY=y1; binY<=y2; binY++){
		histo[chCath]->SetBinContent(binX, binY, efficiency);
		histo[chCath]->SetBinError(binX, binY, efficiencyError);
		histoSlat[chCath]->SetBinContent(binX, binY, efficiencySlat);
		histoSlat[chCath]->SetBinError(binX, binY, efficiencySlatError);
	      }
	    }

	    // Check fired pads (when available)
	    if(padFired[cath]) {
	      histoFired[chCath][localId-1]->Fit("pol0","Q0R","",0., (Float_t)nStrips-1.);
	      fitFunc = histoFired[chCath][localId-1]->GetFunction("pol0");
	      Float_t chi2 = fitFunc->GetChisquare();
	      Float_t ndf = (Float_t)fitFunc->GetNDF();
	      Float_t reducedChi2 = chi2/ndf;
	      if(reducedChi2>kChi2RedMax) {
		isStripOffInBoard[chCath][localId-1] = kTRUE;
		//printf("Chamber = %i Cath = %i Board %i: chi2/NDF = %f\tparam = %f\n", ch, cath, localId, reducedChi2, fitFunc->GetParameter(0));
		fFiredFitHistoList->Add(histoFired[chCath][localId-1]);
	      }
	    }
	  } // loop on cathodes
	} // loop on chambers
      } // loop on boards
    } // loop on crates
  } // loop on DDL

  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    for(Int_t ib=0; ib<kNumOfBoards; ib++){
      fBoardLabelList->Add(boardLabel[chCath][ib]);
    }
  }
  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    for(Int_t ib=0; ib<kNumOfBoards; ib++){
      fBoardLabelList->Add(boardLabelSlat[chCath][ib]);
    }
  }

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
