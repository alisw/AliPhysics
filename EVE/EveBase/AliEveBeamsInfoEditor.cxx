// $Id$
// Author: Stefano Carrazza 2010, CERN, stefano.carrazza@cern.ch

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveBeamsInfo.h"
#include "AliEveBeamsInfoEditor.h"

#include "TVirtualPad.h"
#include "TColor.h"
#include "TGButton.h"
#include "TGColorSelect.h"
#include "TGComboBox.h"
#include "TGDoubleSlider.h"
#include "TGFrame.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGString.h"

//______________________________________________________________________________
// GUI editor for AliEveBeamsInfo.
//

ClassImp(AliEveBeamsInfoEditor)

//______________________________________________________________________________
AliEveBeamsInfoEditor::AliEveBeamsInfoEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fIsMC(0),
  fEventSelection(0),
  fShowEvents(0),
  fSelect(0),
  fButtonPrev(0),
  fButtonNext(0),
  fSetAlpha(0),
  fAlpha(0)
{
  // Constructor.
  MakeTitle("AliEveBeamsInfo");

  // Events selection  
  fShowEvents = new TGCheckButton(this, new TGHotString("&Show information box"));
  fShowEvents->SetState(kButtonDown);
  fShowEvents->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowEventSelection()");
  AddFrame(fShowEvents, new TGLayoutHints(kLHintsExpandX));

  fIsMC = new TGCheckButton(this, new TGHotString("&Data is from simulation (MC)"));
  fIsMC->SetState(kButtonUp);
  fIsMC->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "SwitchDataType()");
  AddFrame(fIsMC, new TGLayoutHints(kLHintsExpandX));

  fEventSelection = new TGGroupFrame(this, "Event filter:", kHorizontalFrame);
  fSelect = new TGComboBox(fEventSelection,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
  fSelect->AddEntry("Show all events",0);
  fSelect->AddEntry("Only Beam 1 events",1);
  fSelect->AddEntry("Only Beam 2 events",2);
  fSelect->AddEntry("Beams 1 & 2 events",3);
  fSelect->Resize(120,22);
  fSelect->Select(0);
  fEventSelection->AddFrame(fSelect, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

  fSelect->Connect("Selected(Int_t)", "AliEveBeamsInfoEditor", this, "SelectEventSelection(Int_t)");

  AddFrame(fEventSelection, new TGLayoutHints(kLHintsExpandX));

  //------

  TGHorizontalFrame *h = new TGHorizontalFrame(this);

  fButtonPrev = new TGTextButton(h, "Previous event");
  h->AddFrame(fButtonPrev, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));
  fButtonPrev->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowPrevEvent()");

  fButtonNext = new TGTextButton(h, "Next event");
  h->AddFrame( fButtonNext, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsExpandX));
  fButtonNext->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowNextEvent()");

  AddFrame(h, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

  fSetAlpha = new TGGroupFrame(this, "Transparency value:", kHorizontalFrame);
  fAlpha = new TGNumberEntry(fSetAlpha,  1.5, 7, -1,
                             TGNumberFormat::kNESRealOne,
                             TGNumberFormat::kNEANonNegative,
                             TGNumberFormat::kNELLimitMinMax,
                             0, 1.5);
  fAlpha->Resize(120,22);
  fSetAlpha->AddFrame(fAlpha, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fAlpha->Connect("ValueSet(Long_t)", "AliEveBeamsInfoEditor", this, "SetAlpha()");

  AddFrame(fSetAlpha, new TGLayoutHints(kLHintsExpandX));
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::SetModel(TObject* obj)
{
  // Set model object.
  fM = dynamic_cast<AliEveBeamsInfo*>(obj);
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::ShowEventSelection()
{
  // Show event selection
  fM->ShowEventSelection(fShowEvents->IsOn());
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::SelectEventSelection(Int_t id)
{
  // Show event selection
  fM->SelectEventSelection(id);
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::ShowPrevEvent()
{
  // Show previous event
  fM->ShowPrevEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::ShowNextEvent()
{
  // Show next event
  fM->ShowNextEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::SwitchDataType()
{
  // Show previous event
  fM->SwitchDataType(fIsMC->IsOn());
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::SetAlpha()
{
  // Set alpha
  fM->SetAlpha(fAlpha->GetNumber());
}


/******************************************************************************/
