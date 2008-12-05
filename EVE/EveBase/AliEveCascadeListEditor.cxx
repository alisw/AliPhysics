// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveCascadeListEditor.h"
#include "AliEveCascade.h"

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
// GUI editor for AliEveCascadeList.
//

ClassImp(AliEveCascadeListEditor)

//______________________________________________________________________________
AliEveCascadeListEditor::AliEveCascadeListEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fMinMaxRCut(0),
  fMinMaxDaughterDCA(0),
  fMinMaxPt(0)
{
  // Constructor.

  MakeTitle("AliEveCascadeList");

  // Create widgets
  // fXYZZ = new TGSomeWidget(this, ...);
  // AddFrame(fXYZZ, new TGLayoutHints(...));
  // fXYZZ->Connect("SignalName()", "Reve::AliEveCascadeListEditor", this, "DoXYZZ()");

   fMinMaxRCut = new TEveGDoubleValuator(this,"Radius:", 130, 0);
   fMinMaxRCut->SetNELength(5);
   fMinMaxRCut->SetLabelWidth(74);
   fMinMaxRCut->Build();
   fMinMaxRCut->GetSlider()->SetWidth(200);
   fMinMaxRCut->SetLimits(0, 100, TGNumberFormat::kNESRealOne);
   fMinMaxRCut->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxRCut()");
   AddFrame(fMinMaxRCut, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxDaughterDCA = new TEveGDoubleValuator(this,"DCA:", 130, 0);
   fMinMaxDaughterDCA->SetNELength(5);
   fMinMaxDaughterDCA->SetLabelWidth(74);
   fMinMaxDaughterDCA->Build();
   fMinMaxDaughterDCA->GetSlider()->SetWidth(200);
   fMinMaxDaughterDCA->SetLimits(0, 1, TGNumberFormat::kNESRealTwo);
   fMinMaxDaughterDCA->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxDaughterDCA()");
   AddFrame(fMinMaxDaughterDCA, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

   fMinMaxPt = new TEveGDoubleValuator(this,"pT:", 130, 0);
   fMinMaxPt->SetNELength(5);
   fMinMaxPt->SetLabelWidth(74);
   fMinMaxPt->Build();
   fMinMaxPt->GetSlider()->SetWidth(200);
   fMinMaxPt->SetLimits(0, 20, TGNumberFormat::kNESRealOne);
   fMinMaxPt->Connect("ValueSet()", "AliEveCascadeListEditor", this, "DoMinMaxPt()");
   AddFrame(fMinMaxPt, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveCascadeListEditor::SetModel(TObject* obj)
{
  // Set model object.

  fM = dynamic_cast<AliEveCascadeList*>(obj);

  // Set values of widgets
  // fXYZZ->SetValue(fM->GetXYZZ());

  fMinMaxRCut->SetValues(fM->fMinRCut, fM->fMaxRCut);
  fMinMaxDaughterDCA->SetValues(fM->fMinDaughterDCA, fM->fMaxDaughterDCA);
  fMinMaxPt->SetValues(fM->fMinPt, fM->fMaxPt);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveCascadeListEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

void AliEveCascadeListEditor::DoMinMaxRCut()
{
  fM->FilterByRadius(fMinMaxRCut->GetMin(), fMinMaxRCut->GetMax());
}

void AliEveCascadeListEditor::DoMinMaxDaughterDCA()
{
  fM->FilterByDaughterDCA(fMinMaxDaughterDCA->GetMin(), fMinMaxDaughterDCA->GetMax());
}

void AliEveCascadeListEditor::DoMinMaxPt()
{
  fM->FilterByPt(fMinMaxPt->GetMin(), fMinMaxPt->GetMax());
}
