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
  fPosCharged(0),
  fNegCharged(0),
  fElectrons(0),
  fMuons(0),
  fPions(0),
  fKaons(0),
  fProtons(0),
  fLabel(0),
  fLabel1(0),
  fThreshold(0),
  fMaxPt(0),
  fSelect(0),
  fParticlesBGAE(0),  
  fTrackSelectionAE(0),
  fPosChargedAE(0),
  fNegChargedAE(0),
  fElectronsAE(0),
  fMuonsAE(0),
  fPionsAE(0),
  fKaonsAE(0),
  fProtonsAE(0),
  fApplyChanges(0),
  fLabelAE(0),
  fLabel1AE(0),
  fThresholdAE(0),
  fMaxPtAE(0),
  fEventControl(0),
  fIsMC(kFALSE),
  fCollisionCandidatesOnly(0)
{
  // Constructor.
  MakeTitle("AliEveLego");

  // Create widgets
  fAllEventsButton = new TGTextButton(this, "Create lego of all events");
  AddFrame(fAllEventsButton, new TGLayoutHints(kLHintsExpandX));
  fAllEventsButton->Connect("Clicked()", "AliEveLegoEditor", this, "DoAllEvents()");

  fParticlesBG = new TGGroupFrame(this, "Particle selection:", kVerticalFrame);
  fPosCharged  = new TGCheckButton(fParticlesBG, new TGHotString("&Positive charged"));
  fNegCharged  = new TGCheckButton(fParticlesBG, new TGHotString("&Negative charged"));
  fElectrons   = new TGCheckButton(fParticlesBG, new TGHotString("&Electrons"));
  fMuons       = new TGCheckButton(fParticlesBG, new TGHotString("&Muons"));
  fPions       = new TGCheckButton(fParticlesBG, new TGHotString("&Pions"));
  fKaons       = new TGCheckButton(fParticlesBG, new TGHotString("&Kaons"));
  fProtons     = new TGCheckButton(fParticlesBG, new TGHotString("&Protons"));

  fPosCharged->SetState(kButtonDown);
  fNegCharged->SetState(kButtonDown);
  fElectrons->SetState(kButtonUp);
  fMuons->SetState(kButtonUp);
  fPions->SetState(kButtonUp);
  fKaons->SetState(kButtonUp);
  fProtons->SetState(kButtonUp);

  fPosCharged->Connect("Clicked()", "AliEveLegoEditor", this, "ShowPosCharge()");
  fNegCharged->Connect("Clicked()", "AliEveLegoEditor", this, "ShowNegCharge()");
  fElectrons->Connect("Clicked()", "AliEveLegoEditor", this, "ShowElectrons()");
  fMuons->Connect("Clicked()", "AliEveLegoEditor", this, "ShowMuons()");
  fPions->Connect("Clicked()", "AliEveLegoEditor", this, "ShowPions()");
  fKaons->Connect("Clicked()", "AliEveLegoEditor", this, "ShowKaons()");
  fProtons->Connect("Clicked()", "AliEveLegoEditor", this, "ShowProtons()");

  fParticlesBG->AddFrame(fPosCharged, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fNegCharged, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fElectrons, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fMuons, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fPions, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fKaons, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->AddFrame(fProtons, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fParticlesBG->SetLayoutManager(new TGVerticalLayout(fParticlesBG));

  AddFrame(fParticlesBG, new TGLayoutHints(kLHintsExpandX));

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


}

/******************************************************************************/

//______________________________________________________________________________
void AliEveLegoEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<AliEveLego*>(obj);
}

//______________________________________________________________________________
void AliEveLegoEditor::DoAllEvents()
{
  fAllEventsButton->SetEnabled(kFALSE);
  CreateAllEventsEditor();
  fM->LoadAllEvents();
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowPosCharge()
{   
  fM->SetParticleType(0);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowNegCharge()
{
  fM->SetParticleType(1);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowElectrons()
{
  fM->SetParticleType(2);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowMuons()
{
  fM->SetParticleType(3);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowPions()
{
  fM->SetParticleType(4);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowKaons()
{
  fM->SetParticleType(5);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowProtons()
{
  fM->SetParticleType(6);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowPosChargeAE()
{
   fM->SetParticleTypeAE(0);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowNegChargeAE()
{
  fM->SetParticleTypeAE(1);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowElectronsAE()
{
  fM->SetParticleTypeAE(2);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowMuonsAE()
{
  fM->SetParticleTypeAE(3);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowPionsAE()
{
  fM->SetParticleTypeAE(4);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowKaonsAE()
{
  fM->SetParticleTypeAE(5);
}

//______________________________________________________________________________
void AliEveLegoEditor::ShowProtonsAE()
{
  fM->SetParticleTypeAE(6);
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
void AliEveLegoEditor::CreateAllEventsEditor()
{
   // create the GUI of all events
   TGVerticalFrame *this2 = this->CreateEditorTabSubFrame("All events style");

   fEventControl = new TGButtonGroup(this2, "Event control:", kVerticalFrame);
   fIsMC = new TGCheckButton(fEventControl, new TGHotString("&Data is from simulation (MC)"));
   fCollisionCandidatesOnly = new TGCheckButton(fEventControl, new TGHotString("&Only collision candidates events"));

   fIsMC->SetState(kButtonUp);
   fCollisionCandidatesOnly->SetState(kButtonUp);

   fIsMC->Connect("Clicked()", "AliEveLegoEditor", this, "DataIsMC()");
   fCollisionCandidatesOnly->Connect("Clicked()", "AliEveLegoEditor", this, "CollisionCandidatesOnly()");
   this2->AddFrame(fEventControl, new TGLayoutHints(kLHintsExpandX));

   fParticlesBGAE = new TGButtonGroup(this2, "Particle selection:", kVerticalFrame);

   fPosChargedAE  = new TGCheckButton(fParticlesBGAE, new TGHotString("&Positive charged"));
   fNegChargedAE  = new TGCheckButton(fParticlesBGAE, new TGHotString("&Negative charged"));
   fElectronsAE   = new TGCheckButton(fParticlesBGAE, new TGHotString("&Electrons"));
   fMuonsAE       = new TGCheckButton(fParticlesBGAE, new TGHotString("&Muons"));
   fPionsAE       = new TGCheckButton(fParticlesBGAE, new TGHotString("&Pions"));
   fKaonsAE       = new TGCheckButton(fParticlesBGAE, new TGHotString("&Kaons"));
   fProtonsAE     = new TGCheckButton(fParticlesBGAE, new TGHotString("&Protons"));

   fPosChargedAE->SetState(kButtonDown);
   fNegChargedAE->SetState(kButtonDown);
   fElectronsAE->SetState(kButtonUp);
   fMuonsAE->SetState(kButtonUp);
   fPionsAE->SetState(kButtonUp);
   fKaonsAE->SetState(kButtonUp);
   fProtonsAE->SetState(kButtonUp);

   fPosChargedAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowPosChargeAE()");
   fNegChargedAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowNegChargeAE()");
   fElectronsAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowElectronsAE()");
   fMuonsAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowMuonsAE()");
   fPionsAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowPionsAE()");
   fKaonsAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowKaonsAE()");
   fProtonsAE->Connect("Clicked()", "AliEveLegoEditor", this, "ShowProtonsAE()");

   this2->AddFrame(fParticlesBGAE, new TGLayoutHints(kLHintsExpandX));

   fApplyChanges = new TGTextButton(this2, "Apply particle selection");
   this2->AddFrame(fApplyChanges, new TGLayoutHints(kLHintsExpandX));
   fApplyChanges->Connect("Clicked()", "AliEveLegoEditor", this, "ApplyChanges()");

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
void AliEveLegoEditor::ApplyChanges()
{
  fM->ApplyParticleTypeSelectionAE();
}

//______________________________________________________________________________
void AliEveLegoEditor::DataIsMC()
{
  fM->SwitchDataType();
}

//______________________________________________________________________________
void AliEveLegoEditor::CollisionCandidatesOnly()
{
  fM->SetCollisionCandidatesOnly();
}
