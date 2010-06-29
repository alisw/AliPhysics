// $Id$
// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveLegoEditor.h"
#include "AliEveLego.h"

#include "TVirtualPad.h"
#include "TColor.h"

#include "TGLabel.h"
#include "TGButton.h"
#include "TGNumberEntry.h"
#include "TGColorSelect.h"
#include "TGDoubleSlider.h"
#include "TGButtonGroup.h"
#include "TGString.h"
#include "TGComboBox.h"
#include "TGFrame.h"

//______________________________________________________________________________
// GUI editor for AliEveLego.
//

ClassImp(AliEveLegoEditor)

//______________________________________________________________________________
AliEveLegoEditor::AliEveLegoEditor(const TGWindow *p, Int_t width, Int_t height,
             UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),

  fAllEventsButton(0),
  fParticlesBG(0),
  fTrackSelection(0),
  fEventSelection(0),
  fRevents(0),
  fLabel(0),
  fLabel1(0),
  fThreshold(0),
  fMaxPt(0),
  fSelect(0),
  fButtonPrev(0),
  fButtonNext(0),

  fParticlesBGAE(0),
  fTrackSelectionAE(0),
  fEventSelectionAE(0),
  fReventsAE(0),
  fLabelAE(0),
  fLabel1AE(0),
  fThresholdAE(0),
  fMaxPtAE(0),
  fSelectAE(0),
  fButtonPrevAE(0),
  fButtonNextAE(0)
{
  // Constructor.
  MakeTitle("AliEveLego");

  // Create widgets
  fAllEventsButton = new TGTextButton(this, "Create lego of all events");
  AddFrame(fAllEventsButton, new TGLayoutHints(kLHintsExpandX));
  fAllEventsButton->Connect("Clicked()", "AliEveLegoEditor", this, "DoAllEvents()");

  fParticlesBG = new TGButtonGroup(this, "Charge selection:", kVerticalFrame);
  fRcharge[0] = new TGRadioButton(fParticlesBG, new TGHotString("&Positive and negative"));
  fRcharge[1] = new TGRadioButton(fParticlesBG, new TGHotString("&Only positive"));
  fRcharge[2] = new TGRadioButton(fParticlesBG, new TGHotString("&Only negative"));
  fRcharge[0]->SetState(kButtonDown);
  AddFrame(fParticlesBG, new TGLayoutHints(kLHintsExpandX));
  fParticlesBG->Connect("Clicked(Int_t)", "AliEveLegoEditor", this, "ShowByCharge(Int_t)");

  fTrackSelection = new TGButtonGroup(this, "Track selection:", kHorizontalFrame);
  fRtracks[0] = new TGRadioButton(fTrackSelection, new TGHotString("&All tracks  "));
  fRtracks[1] = new TGRadioButton(fTrackSelection, new TGHotString("&Primary tracks"));
  fRtracks[0]->SetState(kButtonDown);
  AddFrame(fTrackSelection, new TGLayoutHints(kLHintsExpandX));
  fTrackSelection->Connect("Clicked(Int_t)", "AliEveLegoEditor", this, "ShowByTracks(Int_t)");

  //**************

  TGHorizontalFrame *horz = new TGHorizontalFrame(this);
  AddFrame(horz, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
  fLabel = new TGLabel(horz, "Tracks maximum Pt (GeV): ");
  horz->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

  fMaxPt = new TGNumberEntry(horz, 10000, 7, -1,
                                 TGNumberFormat::kNESRealOne,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax,
                                 0, 10000);

  fMaxPt->Connect("ValueSet(Long_t)", "AliEveLegoEditor", this, "SetMaxPt()");

  horz->AddFrame( fMaxPt, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsCenterY));

  TGHorizontalFrame *horz1 = new TGHorizontalFrame(this);
  AddFrame(horz1, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));
  fLabel1 = new TGLabel(horz1, "Tracks threshold (GeV): ");
  horz1->AddFrame(fLabel1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

  fThreshold = new TGNumberEntry(horz1, 0, 7, -1,
                                 TGNumberFormat::kNESRealOne,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax,
                                 0, 10000);

  fThreshold->Connect("ValueSet(Long_t)", "AliEveLegoEditor", this, "SetThreshold()");
  horz1->AddFrame( fThreshold, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsCenterY));



  //***************

  // Events selection
  fEventSelection = new TGGroupFrame(this, "Event selection:", kHorizontalFrame);

  fRevents = new TGCheckButton(fEventSelection, new TGHotString("&Show info "));
  fRevents->SetState(kButtonUp);
  fRevents->Connect("Clicked()", "AliEveLegoEditor", this, "ShowEventSelection()");
  fEventSelection->AddFrame(fRevents, new TGLayoutHints(kLHintsLeft | kLHintsTop));

  fSelect = new TGComboBox(fEventSelection,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
  fSelect->AddEntry("All events",0);
  fSelect->AddEntry("Beam 1",1);
  fSelect->AddEntry("Beam 2",2);
  fSelect->AddEntry("Beams 1 & 2",3);
  fSelect->Resize(102,22);
  fSelect->Select(0);
  fEventSelection->AddFrame(fSelect, new TGLayoutHints(kLHintsRight | kLHintsExpandX));

  fSelect->Connect("Selected(Int_t)", "AliEveLegoEditor", this, "SelectEventSelection(Int_t)");

  AddFrame(fEventSelection, new TGLayoutHints(kLHintsExpandX));

  //**********

  TGHorizontalFrame *horz3 = new TGHorizontalFrame(this);

  fButtonPrev = new TGTextButton(horz3, "Previous event");
  horz3->AddFrame(fButtonPrev, new TGLayoutHints(kLHintsLeft | kLHintsCenterY | kLHintsExpandX));
  fButtonPrev->Connect("Clicked()", "AliEveLegoEditor", this, "ShowPrevEvent()");

  fButtonNext = new TGTextButton(horz3, "Next event");
  horz3->AddFrame( fButtonNext, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsExpandX));
  fButtonNext->Connect("Clicked()", "AliEveLegoEditor", this, "ShowNextEvent()");

  AddFrame(horz3, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

  //**********


}

/******************************************************************************/

//______________________________________________________________________________
void AliEveLegoEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveLego*>(obj);
}

/******************************************************************************/

// Implements callback/slot methods

//______________________________________________________________________________
// void AliEveLegoEditor::DoXYZZ()
// {
//    // Slot for XYZZ.
//
//    fM->SetXYZZ(fXYZZ->GetValue());
//    Update();
// }

//______________________________________________________________________________
void AliEveLegoEditor::DoAllEvents()
{
  // Slot for XYZZ.
  fAllEventsButton->SetEnabled(kFALSE);
  CreateAllEventsEditor();
  fM->LoadAllEvents();

}

//______________________________________________________________________________
void AliEveLegoEditor::ShowByCharge(Int_t id)
{   
  fM->SetCharge(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowByChargeAE(Int_t id)
{
   fM->SetAllEventsCharge(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::SetMaxPt()
{
   fM->SetMaxPt(fMaxPt->GetNumber());
}

//______________________________________________________________________________
void AliEveLegoEditor::SetMaxPtAE()
{
   fM->SetMaxPtAE(fMaxPtAE->GetNumber());
}

//______________________________________________________________________________
void AliEveLegoEditor::SetThreshold()
{
   fM->SetThreshold(fThreshold->GetNumber());
}

//______________________________________________________________________________
void AliEveLegoEditor::SetThresholdAE()
{
   fM->SetThresholdAE(fThresholdAE->GetNumber());
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowByTracks(Int_t id)
{
   fM->SetTracks(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowByTracksAE(Int_t id)
{
   fM->SetTracksAE(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowByEvents(Int_t id)
{
   fM->SetEvents(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowEventSelection()
{
   fM->SetEventSelection();
}

//______________________________________________________________________________
void AliEveLegoEditor::SelectEventSelection(Int_t id)
{
   fM->SelectEventSelection(id);
}

//______________________________________________________________________________
void AliEveLegoEditor::CreateAllEventsEditor()
{
   // create the GUI of all events
   TGVerticalFrame *this2 = this->CreateEditorTabSubFrame("All events style");

   TGLabel *ftitle = new TGLabel(this2, "AliLego all events ------");
   this2->AddFrame(ftitle, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fParticlesBGAE = new TGButtonGroup(this2, "Charge selection:", kVerticalFrame);
   fRchargeAE[0] = new TGRadioButton(fParticlesBGAE, new TGHotString("&Positive and negative"));
   fRchargeAE[1] = new TGRadioButton(fParticlesBGAE, new TGHotString("&Only positive"));
   fRchargeAE[2] = new TGRadioButton(fParticlesBGAE, new TGHotString("&Only negative"));
   fRchargeAE[0]->SetState(kButtonDown);
   this2->AddFrame(fParticlesBGAE, new TGLayoutHints(kLHintsExpandX));

   fParticlesBGAE->Connect("Clicked(Int_t)", "AliEveLegoEditor", this, "ShowByChargeAE(Int_t)");

   fTrackSelectionAE = new TGButtonGroup(this2, "Track selection:", kHorizontalFrame);
   fRtracksAE[0] = new TGRadioButton(fTrackSelectionAE, new TGHotString("&All tracks  "));
   fRtracksAE[1] = new TGRadioButton(fTrackSelectionAE, new TGHotString("&Primary tracks"));
   fRtracksAE[0]->SetState(kButtonDown);
   this2->AddFrame(fTrackSelectionAE, new TGLayoutHints(kLHintsExpandX));
   fTrackSelectionAE->Connect("Clicked(Int_t)", "AliEveLegoEditor", this, "ShowByTracksAE(Int_t)");

   //**************

   TGHorizontalFrame *horzAE = new TGHorizontalFrame(this2);

   fLabelAE = new TGLabel(horzAE, "Tracks maximum Pt (GeV): ");
   horzAE->AddFrame(fLabelAE, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fMaxPtAE = new TGNumberEntry(horzAE, 10000, 7, -1,
                                  TGNumberFormat::kNESRealOne,
                                  TGNumberFormat::kNEANonNegative,
                                  TGNumberFormat::kNELLimitMinMax,
                                  0, 10000);

   fMaxPtAE->Connect("ValueSet(Long_t)", "AliEveLegoEditor", this, "SetMaxPtAE()");

   horzAE->AddFrame( fMaxPtAE, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsCenterY));
   this2->AddFrame(horzAE, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

   TGHorizontalFrame *horz1AE = new TGHorizontalFrame(this2);

   fLabel1AE = new TGLabel(horz1AE, "Tracks threshold (GeV): ");
   horz1AE->AddFrame(fLabel1AE, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));

   fThresholdAE = new TGNumberEntry(horz1AE, 0, 7, -1,
                                  TGNumberFormat::kNESRealOne,
                                  TGNumberFormat::kNEANonNegative,
                                  TGNumberFormat::kNELLimitMinMax,
                                  0, 10000);

   fThresholdAE->Connect("ValueSet(Long_t)", "AliEveLegoEditor", this, "SetThresholdAE()");
   horz1AE->AddFrame( fThresholdAE, new TGLayoutHints(kLHintsRight | kLHintsNormal | kLHintsCenterY));

   this2->AddFrame(horz1AE, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

}

//______________________________________________________________________________
void AliEveLegoEditor::ShowPrevEvent()
{
   fM->ShowPrevEvent();
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowNextEvent()
{
   fM->ShowNextEvent();
}

