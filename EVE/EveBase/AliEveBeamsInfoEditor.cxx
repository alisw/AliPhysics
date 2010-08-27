// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveBeamsInfoEditor.h"
#include "AliEveBeamsInfo.h"

#include "TVirtualPad.h"
#include "TColor.h"

// Cleanup these includes:
#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include "TGString.h"
#include "TGComboBox.h"
#include "TGFrame.h"


//______________________________________________________________________________
// GUI editor for AliEveBeamsInfo.
//

ClassImp(AliEveBeamsInfoEditor)

//______________________________________________________________________________
AliEveBeamsInfoEditor::AliEveBeamsInfoEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fEventSelection(0),
  fShowEvents(0),
  fSelect(0),
  fButtonPrev(0),
  fButtonNext(0)
{
  // Constructor.
  MakeTitle("AliEveBeamsInfo");

  // Events selection
  fEventSelection = new TGGroupFrame(this, "Event selection:", kHorizontalFrame);

  fShowEvents = new TGCheckButton(fEventSelection, new TGHotString("&Show info "));
  fShowEvents->SetState(kButtonDown);
  fShowEvents->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowEventSelection()");
  fEventSelection->AddFrame(fShowEvents, new TGLayoutHints(kLHintsLeft | kLHintsTop));

  fSelect = new TGComboBox(fEventSelection,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
  fSelect->AddEntry("All events",0);
  fSelect->AddEntry("Beam 1",1);
  fSelect->AddEntry("Beam 2",2);
  fSelect->AddEntry("Beams 1 & 2",3);
  fSelect->Resize(102,22);
  fSelect->Select(0);
  fEventSelection->AddFrame(fSelect, new TGLayoutHints(kLHintsRight | kLHintsExpandX));

  fSelect->Connect("Selected(Int_t)", "AliEveBeamsInfoEditor", this, "SelectEventSelection(Int_t)");

  AddFrame(fEventSelection, new TGLayoutHints(kLHintsExpandX));

  //**********

  TGHorizontalFrame *h = new TGHorizontalFrame(this);

  fButtonPrev = new TGTextButton(h, "Previous event");
  h->AddFrame(fButtonPrev, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));
  fButtonPrev->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowPrevEvent()");

  fButtonNext = new TGTextButton(h, "Next event");
  h->AddFrame( fButtonNext, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsExpandX));
  fButtonNext->Connect("Clicked()", "AliEveBeamsInfoEditor", this, "ShowNextEvent()");

  AddFrame(h, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

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
   fM->ShowEventSelection();
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::SelectEventSelection(Int_t id)
{
   fM->SelectEventSelection(id);
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::ShowPrevEvent()
{
  fM->ShowPrevEvent();
}

//______________________________________________________________________________
void AliEveBeamsInfoEditor::ShowNextEvent()
{
  fM->ShowNextEvent();
}

/******************************************************************************/
