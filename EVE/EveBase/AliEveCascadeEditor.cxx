// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveCascadeEditor.h"
#include "AliEveCascade.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"


//______________________________________________________________________________
// GUI editor for AliEveCascade.
//

ClassImp(AliEveCascadeEditor)

//______________________________________________________________________________
AliEveCascadeEditor::AliEveCascadeEditor(const TGWindow *p, Int_t width, Int_t height,
                               UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fInfoLabel0(0),
  fInfoLabel1(0),
  fXButton(0)
  // Initialize widget pointers to 0
{
  // Constructor.

  MakeTitle("AliEveCascade");

  fInfoLabel0 = new TGLabel(this);
  fInfoLabel0->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel0, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fInfoLabel1 = new TGLabel(this);
  fInfoLabel1->SetTextJustify(kTextLeft);
  AddFrame(fInfoLabel1, new TGLayoutHints(kLHintsLeft|kLHintsExpandX,
                                          8, 0, 2, 0));

  fXButton = new TGTextButton(this, "Detailed View");
  AddFrame(fXButton, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 0, 0));
  fXButton->Connect("Clicked()", "AliEveCascadeEditor", this, "DisplayDetailed()");
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveCascade*>(obj);

  // Set values of widgets
  fInfoLabel0->SetText(Form("Radius = %f, DCA = %f", fM->GetRadius(), fM->GetDaughterDCA()));
  fInfoLabel1->SetText(Form("Pt = %f", fM->GetPt()));
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveCascadeEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

void AliEveCascadeEditor::DisplayDetailed()
{
  printf("Hura!\n");
}
