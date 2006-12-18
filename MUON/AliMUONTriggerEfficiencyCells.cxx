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

#include <fstream>
#include <TRandom.h>
#include "Riostream.h"
#include "TSystem.h"

#include "AliMUONTriggerEfficiencyCells.h"
#include "AliMUONConstants.h"
#include "AliLog.h"

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
  Reset();
}

//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::AliMUONTriggerEfficiencyCells(const char* filename)
:
TObject()
{
///  Constructor using an ASCII file.
  Reset();
  ReadFile(filename);
}


//__________________________________________________________________________
AliMUONTriggerEfficiencyCells::~AliMUONTriggerEfficiencyCells()
{
///  Destructor. Does nothing ;-)
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::GetCellEfficiency(Int_t detElemId, Float_t x, Float_t y, Float_t &eff1, Float_t &eff2)
{
///  Get the efficiencies of the 2 cathodes at a given location (x,y)

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t slat = FindSlatIndex(detElemId);
  TArrayI cell = CellByCoord(detElemId,x,y);
  eff1 = 0.0;
  eff2 = 0.0;
  if(cell.At(0)>=0 && cell.At(1)>=0)
  {
    eff1 = fCellContent[chamber][slat][0][cell.At(0)][cell.At(1)];
    eff2 = fCellContent[chamber][slat][1][cell.At(0)][cell.At(1)];
  }
}


//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::GetCellEfficiency(Int_t detElemId, Int_t localBoard, Float_t &eff1, Float_t &eff2)
{
///  Get the efficiencies of the 2 cathodes at a given local board

  Int_t chamber = FindChamberIndex(detElemId);
  eff1 = fBoardContent[chamber][0][localBoard];
  eff2 = fBoardContent[chamber][1][localBoard];
}


//__________________________________________________________________________
void 
AliMUONTriggerEfficiencyCells::IsTriggered(Int_t detElemId, Float_t x, Float_t y, Bool_t &trig1, Bool_t &trig2)
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
AliMUONTriggerEfficiencyCells::IsTriggered(Int_t detElemId, Int_t localBoard, Bool_t &trig1, Bool_t &trig2)
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
TArrayI AliMUONTriggerEfficiencyCells::CellByCoord(Int_t detElemId, Float_t x, Float_t y)
{
///  Get the efficiencies at a given location.

  Int_t chamber = FindChamberIndex(detElemId);
  Int_t slat = FindSlatIndex(detElemId);
  Int_t cell[2]={-1,-1};
  Float_t maxX = fCellSize[chamber][slat][0]*((Float_t)fCellNumber[chamber][slat][0]);
  Float_t maxY = fCellSize[chamber][slat][1]*((Float_t)fCellNumber[chamber][slat][1]);
  if(x>=0 & x<maxX & y>=0 & y<maxY)
  {
    cell[0] = (Int_t)(x/fCellSize[chamber][slat][0]);
    cell[1] = (Int_t)(y/fCellSize[chamber][slat][1]);
  }
  return TArrayI(2,cell);
}

//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFile(const char* filename)
{
///  Reads a file containing the efficiency map.

  TString fileName = gSystem->ExpandPathName(filename);
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
///  Structure of file containing geometrical efficency
    Int_t datInt=0, detEl=0, chamber=0, rpc=0;
    Float_t datFloat=0.0;
    char dat[50];

    while (file >> dat) {
	    file >> detEl;
	    chamber = FindChamberIndex(detEl);
	    rpc = FindSlatIndex(detEl);
	    file >> dat;
	    for(Int_t i=0; i<2; i++){
        file >> datInt;
        fCellNumber[chamber][rpc][i] = datInt;
        file >> dat;
	    }
	    for(Int_t i=0; i<2; i++){
        file >> datFloat;
        fCellSize[chamber][rpc][i] = datFloat;
        if(i==0)file >> dat;
	    }
	    for(Int_t cath=0; cath<2; cath++){
        file >> dat;
        file >> datInt;
        for(Int_t iy=0; iy<fCellNumber[chamber][rpc][1]; iy++){
          for(Int_t ix=0; ix<fCellNumber[chamber][rpc][0]; ix++){
            file >> datFloat;
            fCellContent[chamber][rpc][cath][ix][iy] = datFloat;
          }
        }
	    }
    }
}

//__________________________________________________________________________
void AliMUONTriggerEfficiencyCells::ReadFileBoards(ifstream &file)
{
///  Structure of file containing local board efficency
    Int_t datInt=0, detEl=0, chamber=0;
    Float_t datFloat=0.0;
    char dat[50];

    while (file >> dat) {
	    file >> detEl;
	    chamber = FindChamberIndex(detEl);
	    //rpc = FindSlatIndex(detEl);
	    for(Int_t cath=0; cath<2; cath++){
		file >> dat;
		file >> datInt;
		for(Int_t board=0; board<fgkNofBoards; board++){
		    file >> datFloat;
		    fBoardContent[chamber][cath][board] = datFloat;
		}
	    }
    }
}

//__________________________________________________________________________
Int_t AliMUONTriggerEfficiencyCells::FindChamberIndex(Int_t detElemId)
{
///  From detElemId to chamber number

  Int_t offset = 100*(AliMUONConstants::NTrackingCh()+1);
  Int_t chamber = (detElemId-offset)/100;
  return chamber;
}


//__________________________________________________________________________
Int_t AliMUONTriggerEfficiencyCells::FindSlatIndex(Int_t detElemId)
{
///  From detElemId to slat index.

  Int_t offset = 100*(AliMUONConstants::NTrackingCh()+1);
  Int_t chamber = FindChamberIndex(detElemId);
  Int_t slat = detElemId-offset-(chamber*100);
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

  for(Int_t chamber=0; chamber<4; chamber++)
  {
    for(Int_t slat=0; slat<18; slat++)
    {
	    for(Int_t cath=0; cath<2; cath++)
      {
        fCellSize[chamber][slat][cath]=0.0;
        fCellNumber[chamber][slat][cath]=0;
        for(Int_t ix=0; ix<fgkNofCells; ix++)
        {
          for(Int_t iy=0; iy<fgkNofCells; iy++)
          {
            fCellContent[chamber][slat][cath][ix][iy]=0.0;
          }
        }
      }
    }
    for(Int_t cath=0; cath<2; cath++)
    {
	for(Int_t board=0; board<fgkNofBoards; board++)
	{
	    fBoardContent[chamber][cath][board]=0.0;
	}
    }
  }
}

