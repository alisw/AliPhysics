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

#include <TBox.h>

#include "AliMUONTriggerGUIboard.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerGUIboard)
/// \endcond

//__________________________________________________________________________
AliMUONTriggerGUIboard::AliMUONTriggerGUIboard(Int_t id, Char_t *name) 
  : TObject(),
    fName(0),
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
    fIsOpen(0)
{
  /// board main constructor

  fName = new TString(name);
  fID   = id;

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
