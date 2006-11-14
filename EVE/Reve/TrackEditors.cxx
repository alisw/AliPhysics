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
#include "TGComboBox.h"

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
  Int_t labelW = 67;

  // --- Limits

  fMaxR = new RGValuator(this, "Max R:", 110, 0);
  fMaxR->SetLabelWidth(labelW);
  fMaxR->SetNELength(6);
  fMaxR->Build();
  fMaxR->SetLimits(0.1, 1000, 100, TGNumberFormat::kNESRealOne);
  fMaxR->SetToolTip("Maximum radius to which the tracks will be drawn.");
  fMaxR->Connect("ValueSet(Double_t)", "Reve::TrackListEditor", this, "DoMaxR()");
  AddFrame(fMaxR, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  fMaxZ = new RGValuator(this, "Max Z:", 110, 0);
  fMaxZ->SetLabelWidth(labelW);
  fMaxZ->SetNELength(6);
  fMaxZ->Build();
  fMaxZ->SetLimits(0.1, 2000, 100, TGNumberFormat::kNESRealOne);
  fMaxZ->SetToolTip("Maximum z-coordinate to which the tracks will be drawn.");
  fMaxZ->Connect("ValueSet(Double_t)", "Reve::TrackListEditor", this, "DoMaxZ()");
  AddFrame(fMaxZ, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

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

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fRnrTracks = new TGCheckButton(f, "Render tracks");
    f->AddFrame(fRnrTracks, new TGLayoutHints(kLHintsLeft, 3, 1, 2, 0));
    fRnrTracks->Connect
      ("Toggled(Bool_t)", "Reve::TrackListEditor", this, "DoRnrTracks()");
    fWidthCombo = new TGLineWidthComboBox(f);
    fWidthCombo->Resize(80, 18);
    f->AddFrame(fWidthCombo, new TGLayoutHints(kLHintsLeft, 8, 1, 0, 0));

    fWidthCombo->Connect
      ("Selected(Int_t)", "Reve::TrackListEditor", this, "DoLineWidth(Int_t)"); 
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 3, 0));
  }

  fRnrMarkers = new TGCheckButton(this, "Render markers");
  AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fRnrMarkers->Connect
    ("Toggled(Bool_t)",
     "Reve::TrackListEditor", this, "DoRnrMarkers()");  

  // --- Kinematics fitting

  fFitDaughters = new TGCheckButton(this, "Fit daughters");
  AddFrame(fFitDaughters, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
  fFitDaughters->Connect("Toggled(Bool_t)","Reve::TrackListEditor", this, "DoFitDaughters()");

  fFitDecay = new TGCheckButton(this, "Fit decay");
  AddFrame(fFitDecay, new TGLayoutHints(kLHintsTop, 3, 1, 2, 0));
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

  fMaxR->SetValue(fTC->GetMaxR());
  fMaxZ->SetValue(fTC->GetMaxZ());
  fMaxOrbits->SetNumber(fTC->GetMaxOrbs());
  fMinAng->SetNumber(fTC->GetMinAng());
  fDelta->SetNumber(fTC->GetDelta());

  fWidthCombo->Select(fTC->GetWidth());

  fRnrTracks->SetState(fTC->GetRnrTracks() ? kButtonDown : kButtonUp);
  fRnrMarkers->SetState(fTC->GetRnrMarkers() ? kButtonDown : kButtonUp);

  fFitDaughters->SetState(fTC->GetFitDaughters() ? kButtonDown : kButtonUp);
  fFitDecay->SetState(fTC->GetFitDecay() ? kButtonDown : kButtonUp);

  fPtRange->SetValues(0.1, 10);
}

/**************************************************************************/

void TrackListEditor::DoMaxR()
{
  fTC->SetMaxR(fMaxR->GetValue());
  Update();
}

void TrackListEditor::DoMaxZ()
{
  fTC->SetMaxZ(fMaxZ->GetValue());
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

void TrackListEditor::DoLineWidth(Int_t width)
{
  fTC->SetWidth(width);
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
  Update();
}
