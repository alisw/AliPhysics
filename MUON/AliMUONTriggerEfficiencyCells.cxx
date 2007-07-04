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
#include "TMatrixF.h"

#include <fstream>
#include <cassert>

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

/// \cond CLASSIMP
ClassImp(AliMUONTriggerEfficiencyCells)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells()
:
TObject()
{
///  Default constructor.
  CheckConstants();
  Reset();
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const char* filename)
:
TObject()
{
///  Constructor using an ASCII file.
  CheckConstants();
  Reset();
  ReadFile(filename);
}


//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::~AliMUONTriggerEfficiencyCells()
{
///  Destructor. Does nothing ;-)
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
  eff1 = fBoardContent[chamber][localBoard-1];
  eff2 = fBoardContent[fgkNchambers+chamber][localBoard-1];
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
    Int_t datInt=0, detEl=0, chamber=0, chCath=0;
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
		for(Int_t board=0; board<AliMpConstants::NofLocalBoards(); board++){
		    file >> datFloat;
		    fBoardContent[chCath][board] = datFloat;
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
    TH1F *histo = 0x0;
    for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    sprintf(histoName, "%sBoardEffChamber%i", cathCode[cath], 11+ch);
	    histo = (TH1F *)file->Get(histoName);
	    if(!(TH1F *)file->Get(histoName)) {
		AliWarning(Form("Can't find histo %s in file %s",histoName, filename));
		continue;
	    }
	    Int_t chCath = fgkNchambers*cath + ch;
	    for(Int_t board=0; board<AliMpConstants::NofLocalBoards(); board++){
		Int_t bin = histo->FindBin(board+1);
		fBoardContent[chCath][board] = histo->GetBinContent(bin);
	    }
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
	fBoardContent[chCath].Set(AliMpConstants::NofLocalBoards());
	fCellSize[chCath].Set(fgkNslats);
	fCellNumber[chCath].Set(fgkNslats);
	for(Int_t slat=0; slat<fgkNslats; slat++){
	    fCellContent[chCath][slat].ResizeTo(fgkNcells,fgkNcells);
	}
    }

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fBoardContent[chCath].Reset();
	fCellSize[chCath].Reset();
	fCellNumber[chCath].Reset();
	for(Int_t slat=0; slat<fgkNslats; slat++){
	    fCellContent[chCath][slat].Zero();
	}
    }
}


