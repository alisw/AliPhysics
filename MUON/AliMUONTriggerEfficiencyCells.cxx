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

#include "AliLog.h"

#include "TRandom.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"

#include "TH2F.h"
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
fNoCountHistoList(0x0)
{
///  Default constructor.
  CheckConstants();
  Reset();
  InitHistos();
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const char* filename)
:
TObject(),
fCountHistoList(0x0),
fNoCountHistoList(0x0)
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
fNoCountHistoList(noCountHistoList)
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
fNoCountHistoList(other.fNoCountHistoList)
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
void AliMUONTriggerEfficiencyCells::ReadFile(const char* filename)
{
///  Reads a file containing the efficiency map.

  TString fileName = gSystem->ExpandPathName(filename);
  if(fileName.EndsWith(".root")){
      ReadHistoBoards(fileName.Data());
      return;
  }

  InitHistos();
  ifstream file(fileName.Data());
  char dat[50];
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
    char dat[50];

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
    char dat[50];

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
void AliMUONTriggerEfficiencyCells::ReadHistoBoards(const char *filename)
{
///  Structure of file (.root) containing local board efficency
    TFile *file = new TFile(filename, "read");
    if(!file) {
	AliWarning(Form("Can't read file %s",filename));
	return;
    }
    char histoName[30];
    char *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

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
///  Sets our internal array contents to zero.

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
void AliMUONTriggerEfficiencyCells::DisplayEfficiency(Bool_t perSlat)
{
  //
  /// Display calculated efficiency.
  //

  const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

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

  Int_t side, col, line, nbx, slat;
  Float_t xCenter, yCenter, zCenter, xWidth, yWidth;
  Float_t x1Label, x2Label, y1Label, y2Label;
  Int_t x1, y1, x2, y2, board=0;
  Char_t name[8], text[200];

  gStyle->SetPalette(1);

  Char_t *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

  Float_t boardsX = 257.00;  // cm
  Float_t boardsY = 307.00;  // cm

  TH2F *histo[fgkNplanes];
  TPaveLabel *boardLabel[fgkNplanes][234];
  assert(kNumOfBoards==234);    
  TArrayI boardsPerColumn[9];
  for(Int_t iLine=0; iLine<9; iLine++){
    boardsPerColumn[iLine].Set(7);
    boardsPerColumn[iLine].Reset();
  }

  Char_t histoName[40], histoTitle[90], labelTxt[5];

  Float_t efficiency, efficiencyError;

  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    for(Int_t ch=0; ch<fgkNchambers; ch++){
      Int_t chCath = fgkNchambers*cath + ch;
      sprintf(histoName, "%sChamber%i", cathCode[cath], 11+ch);
      sprintf(histoTitle, "Chamber %i: efficiency %s", 11+ch, cathCode[cath]);
      histo[chCath] = new TH2F(histoName, histoTitle, (Int_t)boardsX, -boardsX, boardsX, (Int_t)boardsY, -boardsY, boardsY);
      histo[chCath]->SetXTitle("X (cm)");
      histo[chCath]->SetYTitle("Y (cm)");
    }
  }

  TString mapspath = gSystem->Getenv("ALICE_ROOT");
  mapspath.Append("/MUON/data");

  sprintf(text,"%s/guimapp11.txt",mapspath.Data());

  FileStat_t fs;
  if(gSystem->GetPathInfo(text,fs)){
    AliWarning(Form("Map file %s not found. Nothing done",text));
    return;
  }

  FILE *fmap = fopen(text,"r");

  for (Int_t ib = 0; ib < kNumOfBoards; ib++) {
    fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);
    if(side==0) continue;
    boardsPerColumn[line-1][col-1]++;
  }

  rewind(fmap);

  for (Int_t ib = 0; ib < kNumOfBoards; ib++) {
    fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);

    board=0;
    for(Int_t iCol=1; iCol<=col; iCol++){
      Int_t lastLine = 9;
      if(iCol==col) lastLine = line-1;
      for(Int_t iLine=1; iLine<=lastLine; iLine++){
	board += boardsPerColumn[iLine-1][iCol-1];
      }
    }
    if(side==0) board += kNumOfBoards/2;
    board += nbx - 1;

    slat = (line+13)%fgkNslats;
    if(side==0) slat = 14-line;

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
      x1 = histo[chCath]->GetXaxis()->FindBin(xCenter-xWidth/2.)+1;
      y1 = histo[chCath]->GetYaxis()->FindBin(yCenter-yWidth/2.)+1;
      x2 = histo[chCath]->GetXaxis()->FindBin(xCenter+xWidth/2.)-1;
      y2 = histo[chCath]->GetYaxis()->FindBin(yCenter+yWidth/2.)-1;

      x1Label = xCenter-xWidth/2.;
      y1Label = yCenter-yWidth/2.;
      x2Label = xCenter+xWidth/2.;
      y2Label =  yCenter+yWidth/2.;
      sprintf(labelTxt,"%3d", board+1);
      if(perSlat){
	x1Label = 140.;
	x2Label = x1Label + 40.;
	y1Label = -285. + ((Float_t)(line - 1)) * 68;
	y2Label = y1Label + 34.;
	if(side==0){
	  x1Label = -x2Label;
	  x2Label = x1Label + 40.;
	}
	sprintf(labelTxt,"%2d", slat);
      }
    
      boardLabel[chCath][board] = new TPaveLabel(x1Label, y1Label, x2Label, y2Label, labelTxt);
      boardLabel[chCath][board]->SetFillStyle(0);
      boardLabel[chCath][board]->SetBorderSize(0);

      Int_t histoBin = board+1;
      if(perSlat) histoBin = slat+1;

      if(!perSlat){
	efficiency = fBoardEfficiency[chCath]->GetBinContent(histoBin);
	efficiencyError = fBoardEfficiency[chCath]->GetBinError(histoBin);
      }
      else {
	efficiency = fSlatEfficiency[chCath]->GetBinContent(histoBin);
	efficiencyError = fSlatEfficiency[chCath]->GetBinError(histoBin);
      }

      for(Int_t binX=x1; binX<=x2; binX++){
	for(Int_t binY=y1; binY<=y2; binY++){
	  histo[chCath]->SetBinContent(binX, binY, efficiency);
	  histo[chCath]->SetBinError(binX, binY, efficiencyError);
	}
      }
    }
  }

  fclose(fmap);

  TCanvas *can[8];
  for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
    sprintf(histoName, "%sCan", histo[chCath]->GetName());
    sprintf(histoTitle, "%s", histo[chCath]->GetTitle());
    can[chCath] = new TCanvas(histoName, histoTitle, 100+10*chCath, 10*chCath, 700, 700);
    can[chCath]->SetRightMargin(0.14);
    can[chCath]->SetLeftMargin(0.12);
    histo[chCath]->GetZaxis()->SetRangeUser(0.,1.);
    histo[chCath]->GetYaxis()->SetTitleOffset(1.4);
    histo[chCath]->SetStats(kFALSE);
    histo[chCath]->Draw("COLZ");
    for (Int_t board = 0; board < kNumOfBoards; board++) {
      boardLabel[chCath][board]->Draw("same");
    }
  }
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
  return kTRUE;
}
