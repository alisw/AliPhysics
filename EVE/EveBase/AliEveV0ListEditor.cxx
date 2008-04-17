// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveV0ListEditor.h"
#include "AliEveV0.h"

#include "TEveGValuators.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"


//______________________________________________________________________________
// GUI editor for AliEveV0List.
//

ClassImp(AliEveV0ListEditor)

//______________________________________________________________________________
AliEveV0ListEditor::AliEveV0ListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fMinMaxRCut(0)
{
  // Constructor.

  MakeTitle("AliEveV0List");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::AliEveV0ListEditor", this, "DoXYZZ()");

   fMinMaxRCut = new TEveGDoubleValuator(this,"Radius:", 130, 0);
   fMinMaxRCut->SetNELength(5);
   fMinMaxRCut->SetLabelWidth(74);
   fMinMaxRCut->Build();
   fMinMaxRCut->GetSlider()->SetWidth(200);
   fMinMaxRCut->SetLimits(0, 100, TGNumberFormat::kNESRealOne);
   fMinMaxRCut->Connect("ValueSet()", "AliEveV0ListEditor", this, "DoMinMaxRCut()");
   AddFrame(fMinMaxRCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveV0ListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveV0List*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxRCut->SetValues(fM->fMinRCut, fM->fMaxRCut);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveV0ListEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

void AliEveV0ListEditor::DoMinMaxRCut()
{
  fM->FilterByRadius(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}
