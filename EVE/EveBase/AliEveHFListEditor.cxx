// @(#)root/eve:$Id$
// Main author: Davide Caffarri 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHFListEditor.h"
#include "AliEveHF.h"

#include "TEveGValuators.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include "TGComboBox.h"
#include "TGLabel.h"
#include "TGSlider.h"

//______________________________________________________________________________
// GUI editor for AliEveHFList.
//

ClassImp(AliEveHFListEditor)

//______________________________________________________________________________
AliEveHFListEditor::AliEveHFListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fMinMaxPt(0),
  fMinMaxCosPointingAngle(0),
  fValueInvMass(0)

{
  // Constructor.

  MakeTitle("AliEveHFList");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::AliEveHFListEditor", this, "DoXYZZ()");

   fMinMaxPt = new TEveGDoubleValuator(this,"pT:", 80, 0);
   fMinMaxPt->SetNELength(5);
   fMinMaxPt->SetLabelWidth(120);
   fMinMaxPt->Build();
   fMinMaxPt->GetSlider()->SetWidth(200);
   fMinMaxPt->SetLimits(0, 20, TGNumberFormat::kNESRealOne);
   fMinMaxPt->Connect("ValueSet()", "AliEveHFListEditor", this, "DoMinMaxPt()");
   AddFrame(fMinMaxPt, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxCosPointingAngle = new TEveGDoubleValuator(this,"Cos Pointing Angle:", 80, 0);
   fMinMaxCosPointingAngle->SetNELength(5);
   fMinMaxCosPointingAngle->SetLabelWidth(120);
   fMinMaxCosPointingAngle->Build();
   fMinMaxCosPointingAngle->GetSlider()->SetWidth(200);
   fMinMaxCosPointingAngle->SetLimits(0, 1, TGNumberFormat::kNESRealOne);
   fMinMaxCosPointingAngle->Connect("ValueSet()", "AliEveHFListEditor", this, "DoMinMaxCosPointingAngle()");
   AddFrame(fMinMaxCosPointingAngle, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));


   fValueInvMass = new TEveGDoubleValuator(this, "Delta Invariant Mass:", 80, 0);
   fValueInvMass->SetNELength(5);
   fValueInvMass->SetLabelWidth(120);
   fValueInvMass->Build();
   fValueInvMass->GetSlider()->SetWidth(200);
   fValueInvMass->SetLimits(0,1, TGNumberFormat::kNESRealTwo);
   fValueInvMass->Connect("ValueSet()", "AliEveHFListEditor", this, "DoMinMaxInvMass()");
   AddFrame(fValueInvMass, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));


}



/******************************************************************************/

//______________________________________________________________________________
void AliEveHFListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = static_cast<AliEveHFList*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxPt->SetValues(fM->fMinPt, fM->fMaxPt);
  fMinMaxCosPointingAngle->SetValues(fM->fMinCosPointingAngle, fM->fMaxCosPointingAngle);
  fValueInvMass->SetValues(0, fM->fDeltaInvariantMass);
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



void AliEveHFListEditor::DoMinMaxPt()
{
  fM->FilterByPt(fMinMaxPt->GetMin(), fMinMaxPt->GetMax());
}

void AliEveHFListEditor::DoMinMaxCosPointingAngle()
{
  fM->FilterByCosPointingAngle(fMinMaxCosPointingAngle->GetMin(), fMinMaxCosPointingAngle->GetMax());
}

void AliEveHFListEditor::DoMinMaxInvMass()
{
  fM->FilterByInvariantMass(0, fValueInvMass->GetMax());
}
