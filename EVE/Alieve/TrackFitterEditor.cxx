// $Header$

#include "TrackFitterEditor.h"
#include <Alieve/TrackFitter.h>

#include <TGButton.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TrackFitterEditor
//

ClassImp(TrackFitterEditor)

TrackFitterEditor::TrackFitterEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fFit(0),
  fReset(0),
  fStart(0),
  fStop(0),
  fGraph(0)
{
  // Constructor.

  MakeTitle("TrackFitter");

  fStart = new TGTextButton(this, "Start");
  AddFrame(fStart, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 1, 1, 1));
  fStart->Connect("Clicked()",
                  "Alieve::TrackFitterEditor", this, "DoStart()");

  fFit = new TGTextButton(this, "Fit");
  AddFrame(fFit, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 1, 1, 1));
  fFit->Connect("Clicked()",
                "Alieve::TrackFitterEditor", this, "DoFit()");

  fReset = new TGTextButton(this, "Reset");
  AddFrame(fReset, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 1, 1, 1));
  fReset->Connect("Clicked()",
                  "Alieve::TrackFitterEditor", this, "DoReset()");

  fStop = new TGTextButton(this, "Stop");
  AddFrame(fStop, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 4, 1, 1, 1));
  fStop->Connect("Clicked()",
                 "Alieve::TrackFitterEditor", this, "DoStop()");

  fGraph = new TGTextButton(this, " RiemanGraph ");
  AddFrame(fGraph, new TGLayoutHints(kLHintsLeft, 4, 2, 4, 1));
  fGraph->Connect("Clicked()",
                 "Alieve::TrackFitterEditor", this, "DoGraph()");
 }

/**************************************************************************/

void TrackFitterEditor::SetModel(TObject* obj)
{ 
  // Set model object.

  fM = dynamic_cast<TrackFitter*>(obj);

  if(fM->GetConnected())
  {
    fStart->SetState(kButtonDisabled);
    fStop->SetState(kButtonUp);
  }
  else
  { 
    fStop->SetState(kButtonDisabled);
    fStart->SetState(kButtonEngaged);
    fStart->SetState(kButtonUp);
  }
}

void TrackFitterEditor::DoFit()
{
  // Fit slot.

  fM->FitTrack();
  Update();
}

void TrackFitterEditor::DoReset()
{
  // Reset slot.

  fM->Reset();
  Update();
}

void TrackFitterEditor::DoStart()
{
  // Start selection slot.

  fM->Start();
  fStart->SetState(kButtonDisabled);
  fStop->SetState(kButtonUp);
}

void TrackFitterEditor::DoStop()
{
  // Stop selection slot.

  fM->Stop();
  fStop->SetState(kButtonDisabled);
  fStart->SetState(kButtonUp);
}

/**************************************************************************/

void TrackFitterEditor::DoGraph()
{
  // Draw graph slot.

  fM->DrawRiemanGraph();
  Update();
}
