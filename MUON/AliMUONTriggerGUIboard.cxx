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

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerGUIboard
///
/// Single trigger board object with geometry information, strips and digits
///
/// \author Bogdan Vulpescu, LPC Clermont-Ferrand
//-----------------------------------------------------------------------------

#include <TClonesArray.h>
#include <TBox.h>
#include <TMath.h>

#include "AliMUONGeometryTransformer.h"

#include "AliMUONTriggerGUIboard.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerGUIboard)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerGUIboard::AliMUONTriggerGUIboard() 
  : TObject(),
    fName(0),
    fCrateName(0),
    fID(-1),
    fStatus(0),
    fPosition(0),
    fYOver(0),
    fXSix(0),
    fXSiy1(0),
    fXSiy2(0),
    fYSix1(0),
    fYSix2(0),
    fYSiy(0),
    fDetElemId(0),
    fIdCircuit(-1),
    fIsOpen(0),
    fNPadsX(),
    fNPadsY(),
    fPadsX(),
    fPadsY()
{
  /// board main constructor

  fName = new TString("");
  fCrateName = new TString("");

  for (Int_t i = 0; i < kNMT; i++) {
    fXCenter[i] = 0.;
    fYCenter[i] = 0.;
    fZCenter[i] = 0.;
    fXWidth[i]  = 0.;
    fYWidth[i]  = 0.;
    for (Int_t is = 0; is < kNS; is++) {
      fXDig[i][is] = 0;
      fYDig[i][is] = 0;
      fXDigBox[i][is] = new TBox(0,0,0,0);
      fYDigBox[i][is] = new TBox(0,0,0,0);
      fXDigBox[i][is]->SetBit(kCannotPick);
      fYDigBox[i][is]->SetBit(kCannotPick);
      fXDigBox[i][is]->SetFillStyle(1001);
      fYDigBox[i][is]->SetFillStyle(1001);
      fXDigBox[i][is]->SetFillColor(4);
      fYDigBox[i][is]->SetFillColor(4);
    }
  }

  fXSix  = -1;
  fXSiy1 = -1;
  fXSiy2 = -1;

  fYSix1 = -1;
  fYSix2 = -1;
  fYSiy  = -1;

  fDetElemId = -1;
  fIdCircuit = -1;

  fIsOpen = kFALSE;

  fYOver    = 0;
  fPosition = 0;

  for (Int_t i = 0; i < kNMT; i++) {
    fPadsX[i] = new TClonesArray("AliMpPad",16); fNPadsX[i] = 0;
    fPadsY[i] = new TClonesArray("AliMpPad",16); fNPadsY[i] = 0;
  }

}

//__________________________________________________________________________
AliMUONTriggerGUIboard::~AliMUONTriggerGUIboard() 
{
  /// board destructor

  delete fName;

  for (Int_t imt = 0; imt < kNMT; imt++) {
    for (Int_t is = 0; is < kNS; is++) {
      delete fXDigBox[imt][is];
      delete fYDigBox[imt][is];
    }
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::SetXDigBox(Int_t imt, Int_t is, Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  /// set coordinates of "is" x-strip box in chamber "imt"

  fXDigBox[imt][is]->SetX1(x1);
  fXDigBox[imt][is]->SetY1(y1);
  fXDigBox[imt][is]->SetX2(x2);
  fXDigBox[imt][is]->SetY2(y2);

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::SetYDigBox(Int_t imt, Int_t is, Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  /// set coordinates of "is" y-strip box in chamber "imt"

  fYDigBox[imt][is]->SetX1(x1);
  fYDigBox[imt][is]->SetY1(y1);
  fYDigBox[imt][is]->SetX2(x2);
  fYDigBox[imt][is]->SetY2(y2);

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::ClearXDigits()
{
  /// delete the set x-digits

  for (Int_t imt = 0; imt < kNMT; imt++) {
    for (Int_t is = 0; is < kNS; is++) {
      fXDig[imt][is] = 0;
    }
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::ClearYDigits()
{
  /// delete the set y-digits

  for (Int_t imt = 0; imt < kNMT; imt++) {
    for (Int_t is = 0; is < kNS; is++) {
      fYDig[imt][is] = 0;
    }
  }

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::MakeGeometry()
{
  /// create the display geometry from the mapping pads

  AliMpPad *pad;

  // circuit number and manu channel (from x-strips)
  for (Int_t ich = 0; ich < kNMT; ich++) {
    if (fNPadsX[ich]) {
      pad = (AliMpPad*)fPadsX[ich]->At(0);
      fIdCircuit = pad->GetLocalBoardId(0);
      break;
    }
  }

  // position index
  if (fName->Length()) {
    if (fName->Contains("12")) fPosition = 1;
    if (fName->Contains("34")) fPosition = 2;
    if (fName->Contains("56")) fPosition = 3;
    if (fName->Contains("78")) fPosition = 4;
  }

  // position index for common y-strip boards
  for (Int_t ich = 0; ich < kNMT; ich++) {
    if (fNPadsY[ich]) {
      pad = (AliMpPad*)fPadsY[ich]->At(0);
      fYOver = pad->GetNofLocations();
      break;
    }
  }

  // pad indices
  Int_t padxIx = -1, padxIy1 = +999, padxIy2 = -999;
  Int_t padyIy = -1, padyIx1 = +999, padyIx2 = -999;
  for (Int_t ip = 0; ip < fNPadsX[0]; ip++) {
    pad = (AliMpPad*)fPadsX[0]->At(ip);
    padxIx = pad->GetIx();
    padxIy1 = TMath::Min(padxIy1,pad->GetIy());
    padxIy2 = TMath::Max(padxIy2,pad->GetIy());
  }
  for (Int_t ip = 0; ip < fNPadsY[0]; ip++) {
    pad = (AliMpPad*)fPadsY[0]->At(ip);
    padyIy = pad->GetIy();
    padyIx1 = TMath::Min(padyIx1,pad->GetIx());
    padyIx2 = TMath::Max(padyIx2,pad->GetIx());
  }
  fXSix  = padxIx;
  fXSiy1 = padxIy1;
  fXSiy2 = padxIy2;
  fYSiy  = padyIy;
  fYSix1 = padyIx1;
  fYSix2 = padyIx2;

  // position and dimension

  AliMUONGeometryTransformer transformer;
  transformer.LoadGeometryData("transform.dat");

  Float_t minX, maxX, minY, maxY;
  Float_t dx, dy;
  Float_t xloc, yloc, xglo, yglo, zglo;
  for (Int_t ich = 0; ich < kNMT; ich++) {
    minX = +9999; maxX = -9999;
    minY = +9999; maxY = -9999;
    for (Int_t ix = 0; ix < fNPadsX[ich]; ix++) {
      pad = (AliMpPad*)fPadsX[ich]->At(ix);
      xloc = pad->GetPositionX();
      yloc = pad->GetPositionY();
      dx = pad->GetDimensionX();
      dy = pad->GetDimensionY();
      transformer.Local2Global((11+ich)*100+GetDetElemId(), xloc, yloc, 0, xglo, yglo, zglo);
      minX = TMath::Min(minX,(xglo-dx));
      maxX = TMath::Max(maxX,(xglo+dx));
      minY = TMath::Min(minY,(yglo-dy));
      maxY = TMath::Max(maxY,(yglo+dy));
    }
    fXCenter[ich] = 0.5*(minX+maxX);
    fYCenter[ich] = 0.5*(minY+maxY);
    fZCenter[ich] = zglo;
    fXWidth[ich]  = maxX-minX;
    fYWidth[ich]  = maxY-minY;
    // truncate to same precision as in the old guimap files
    fXCenter[ich] = 0.01*TMath::Nint(fXCenter[ich]*100.0);
    fYCenter[ich] = 0.01*TMath::Nint(fYCenter[ich]*100.0);
    fXWidth[ich] = 0.01*TMath::Nint(fXWidth[ich]*100.0);
    fYWidth[ich] = 0.01*TMath::Nint(fYWidth[ich]*100.0);

  }

  // delete the pads arrays
  for (Int_t ich = 0; ich < kNMT; ich++) {
    delete fPadsX[ich]; fNPadsX[ich] = 0;
    delete fPadsY[ich]; fNPadsY[ich] = 0;
  }
  
}

//__________________________________________________________________________
Int_t AliMUONTriggerGUIboard::GetLine() const
{
  /// get detector side
  if (fName->Length() >= 5) {
    const Char_t *name = fName->Data();
    TString sline = TString(name[4]);
    return sline.Atoi();
  }

  return -1;

}

//__________________________________________________________________________
Int_t AliMUONTriggerGUIboard::GetCol() const
{
  /// get detector side
  if (fName->Length() >= 5) {
    const Char_t *name = fName->Data();
    TString scol = TString(name[2]);
    return scol.Atoi();
  }

  return -1;

}

//__________________________________________________________________________
Int_t AliMUONTriggerGUIboard::GetSide() const
{
  /// get detector side
  if (fName->Length() >= 5) {
    const Char_t *name = fName->Data();
    if (!strncmp(name,"L",1)) return 0;
    if (!strncmp(name,"R",1)) return 1;
  }

  return -1;

}

//__________________________________________________________________________
void AliMUONTriggerGUIboard::PrintBoard() const
{
  /// print information on this board

  printf("Name: %s Id %3d Circ %3d DetElemId %2d Pos %1d YOver %1d\n",GetBoardName(),GetNumber(),GetIdCircuit(),GetDetElemId(),GetPosition(),GetYOver());
  printf("NStrips: X %2d Y %2d \n",GetNStripX(),GetNStripY());
  printf("Pad indices: X: ix %3d iy1 %3d iy2 %3d \n",GetXSix(),GetXSiy1(),GetXSiy2());
  printf("Pad indices: Y: iy %3d ix1 %3d ix2 %3d \n",GetYSiy(),GetYSix1(),GetYSix2());
  printf("Position and dimension:\n");
  for (Int_t imt = 0; imt < 4; imt++) {
    printf("MT=%1d: X %9.4f Y %9.4f Z %10.4f \n",imt,GetXCenter(imt),GetYCenter(imt),GetZCenter(imt));
    printf("      DX %7.4f DY %7.4f \n",GetXWidth(imt),GetYWidth(imt));
  }

}

