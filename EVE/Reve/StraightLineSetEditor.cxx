// $Header$

#include "StraightLineSetEditor.h"
#include <Reve/StraightLineSet.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// StraightLineSetEditor
//

ClassImp(StraightLineSetEditor)

StraightLineSetEditor::StraightLineSetEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0)
  // Initialize widget pointers to 0
{
  MakeTitle("StraightLineSet");

  TGHorizontalFrame* frame = new TGHorizontalFrame(this);

  fRnrMarkers = new TGCheckButton(frame, "RnrMarkers");
  frame->AddFrame(fRnrMarkers, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));
  fRnrMarkers->Connect
    ("Toggled(Bool_t)",
     "Reve::StraightLineSetEditor", this, "DoRnrMarkers()");

  fRnrLines = new TGCheckButton(frame, "RnrLines");
  frame->AddFrame(fRnrLines, new TGLayoutHints(kLHintsLeft, 2, 1, 1, 1));
  fRnrLines->Connect
    ("Toggled(Bool_t)",
     "Reve::StraightLineSetEditor", this, "DoRnrLines()");

  AddFrame(frame, new TGLayoutHints(kLHintsTop, 0, 0, 0, 0));
}

StraightLineSetEditor::~StraightLineSetEditor()
{}

/**************************************************************************/

void StraightLineSetEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<StraightLineSet*>(obj);

  // Set values of widgets
  fRnrMarkers->SetState(fM->GetRnrMarkers() ? kButtonDown : kButtonUp);
  fRnrLines->SetState(fM->GetRnrLines() ? kButtonDown : kButtonUp);
}

/**************************************************************************/

// Implements callback/slot methods

void StraightLineSetEditor::DoRnrMarkers()
{
  fM->SetRnrMarkers(fRnrMarkers->IsOn());
   Update();
}

void StraightLineSetEditor::DoRnrLines()
{
  fM->SetRnrLines(fRnrLines->IsOn());
  Update();
}
