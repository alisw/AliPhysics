// $Header$

#include "TrackEditors.h"
#include <Reve/Track.h>

#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// TrackListEditor
//

ClassImp(TrackListEditor)

TrackListEditor::TrackListEditor(const TGWindow *p,
                                 Int_t width, Int_t height,
                                 UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),

  fTC (0),

  fMaxR(0),
  fMaxZ(0),
  fMaxOrbits(0),
  fMinAng(0),
  fDelta(0),

  fRnrTracks(0),
  fRnrMarkers(0),

  fFitDaughters(0),
  fFitDecay(0),

  fPtRange(0)
{
  MakeTitle("TrackList");

    // --- Limits

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Max R:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 25, 2, 1, 1));
    fMaxR = new TGNumberEntry(f, 0., 6, -1, 
			      TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
			      TGNumberFormat::kNELLimitMinMax, 0.1, 2000.0);
    fMaxR->GetNumberEntry()->SetToolTipText("Maximum radius to which the tracks will be drawn.");
    f->AddFrame(fMaxR, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMaxR->Associate(f);
    fMaxR->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMaxR()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Max Z:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 24, 2, 1, 1));
    fMaxZ = new TGNumberEntry(f, 0., 6, -1, 
			      TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
			      TGNumberFormat::kNELLimitMinMax, 0.1, 2000.0);
    fMaxZ->GetNumberEntry()->SetToolTipText("Maximum z-coordinate to which the tracks will be drawn.");
    f->AddFrame(fMaxZ, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMaxZ->Associate(f);
    fMaxZ->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMaxZ()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Max Orbits:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 0, 2, 1, 1));
    fMaxOrbits = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				   TGNumberFormat::kNELLimitMinMax, 0.1, 100.0);
    fMaxOrbits->GetNumberEntry()->SetToolTipText("Maximal angular path of tracks' orbits (1 ~ 2Pi).");
    f->AddFrame(fMaxOrbits, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMaxOrbits->Associate(f);
    fMaxOrbits->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMaxOrbits()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Min Angle:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 3, 2, 1, 1));
    fMinAng = new TGNumberEntry(f, 0., 6, -1, 
				TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELLimitMinMax, 1, 180.0);
    fMinAng->GetNumberEntry()->SetToolTipText("Minimal angular step between two helix points.");
    f->AddFrame(fMinAng, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMinAng->Associate(f);
    fMinAng->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoMinAng()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Delta:");
    f->AddFrame(l, new TGLayoutHints(kLHintsTop | kLHintsCenterY, 32, 2, 1, 1));
    fDelta = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				   TGNumberFormat::kNELLimitMinMax, 0.001, 100.0);
    fDelta->GetNumberEntry()->SetToolTipText("Maximal error at the mid-point of the line connecting to helix points.");
    f->AddFrame(fDelta, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fDelta->Associate(f);
    fDelta->Connect("ValueSet(Long_t)", "Reve::TrackListEditor", this, "DoDelta()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  // --- Rendering control

  fRnrTracks = new TGCheckButton(this, "Render tracks");
  AddFrame(fRnrTracks, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrTracks->Connect
    ("Toggled(Bool_t)",
     "Reve::TrackListEditor", this, "DoRnrTracks()");

  fRnrMarkers = new TGCheckButton(this, "Render markers");
  AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrMarkers->Connect
    ("Toggled(Bool_t)",
     "Reve::TrackListEditor", this, "DoRnrMarkers()");  

  // --- Kinematics fitting

  fFitDaughters = new TGCheckButton(this, "Fit daughters");
  AddFrame(fFitDaughters, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fFitDaughters->Connect("Toggled(Bool_t)","Reve::TrackListEditor", this, "DoFitDaughters()");

  fFitDecay = new TGCheckButton(this, "Fit decay");
  AddFrame(fFitDecay, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fFitDecay->Connect("Toggled(Bool_t)","Reve::TrackListEditor", this, "DoFitDecay()");  

  // --- Selectors

  fPtRange = new RGDoubleValuator(this,"Pt Range", 200, 0);
  fPtRange->SetNELength(6);
  fPtRange->Build();
  fPtRange->GetSlider()->SetWidth(224);
  fPtRange->SetLimits(0.1, 10, TGNumberFormat::kNESRealTwo);
  fPtRange->Connect("ValueSet()",
                    "Reve::TrackListEditor", this, "DoPtRange()");
  AddFrame(fPtRange, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
}

TrackListEditor::~TrackListEditor()
{}

/**************************************************************************/

void TrackListEditor::SetModel(TObject* obj)
{
  fTC = dynamic_cast<TrackList*>(obj);

  fMaxR->SetNumber(fTC->GetMaxR());
  fMaxZ->SetNumber(fTC->GetMaxZ());
  fMaxOrbits->SetNumber(fTC->GetMaxOrbs());
  fMinAng->SetNumber(fTC->GetMinAng());
  fDelta->SetNumber(fTC->GetDelta());

  fRnrTracks->SetState(fTC->GetRnrTracks() ? kButtonDown : kButtonUp);
  fRnrMarkers->SetState(fTC->GetRnrMarkers() ? kButtonDown : kButtonUp);

  fFitDaughters->SetState(fTC->GetFitDaughters() ? kButtonDown : kButtonUp);
  fFitDecay->SetState(fTC->GetFitDecay() ? kButtonDown : kButtonUp);

  fPtRange->SetValues(0.1, 10);
}

/**************************************************************************/

void TrackListEditor::DoMaxR()
{
  Double_t maxr = fMaxR->GetNumber();
  fTC->SetMaxR(maxr);
  Update();
}

void TrackListEditor::DoMaxZ()
{
  fTC->SetMaxZ(fMaxZ->GetNumber());
  Update();
}

void TrackListEditor::DoMaxOrbits()
{
  fTC->SetMaxOrbs(fMaxOrbits->GetNumber());
  Update();
}

void TrackListEditor::DoMinAng()
{
  fTC->SetMinAng(fMinAng->GetNumber());
  Update();
}

void TrackListEditor::DoDelta()
{
  fTC->SetDelta(fDelta->GetNumber());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoRnrTracks()
{
  fTC->SetRnrTracks(fRnrTracks->IsOn());
  Update();
}

void TrackListEditor::DoRnrMarkers()
{
  fTC->SetRnrMarkers(fRnrMarkers->IsOn());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoFitDaughters()
{
  fTC->SetFitDaughters(fFitDaughters->IsOn());
  Update();
}

void TrackListEditor::DoFitDecay()
{
  fTC->SetFitDecay(fFitDecay->IsOn());
  Update();
}

/**************************************************************************/

void TrackListEditor::DoPtRange()
{
  fTC->SelectByPt(fPtRange->GetMin(), fPtRange->GetMax());
}
