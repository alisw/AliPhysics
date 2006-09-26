// $Header$

#include "TPCSector2DEditor.h"
#include <Alieve/TPCSector2D.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCSector2DEditor
//

ClassImp(TPCSector2DEditor)

TPCSector2DEditor::TPCSector2DEditor(const TGWindow *p,
                                     Int_t width, Int_t height,
				     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fShowMax(0), fAverage(0), fUseTexture(0)
{
  MakeTitle("TPCSector2D");

  fUseTexture = new TGCheckButton(this, "UseTexture");
  AddFrame(fUseTexture, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fUseTexture->Connect
    ("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoUseTexture()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fShowMax = new TGCheckButton(f, "ShowMax");
    f->AddFrame(fShowMax, new TGLayoutHints(kLHintsLeft, 3, 16, 1, 0));
    fShowMax->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoShowMax()");
    fAverage = new TGCheckButton(f, "Average");
    f->AddFrame(fAverage, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 0));
    fAverage->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoAverage()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
}

TPCSector2DEditor::~TPCSector2DEditor()
{}

/**************************************************************************/

void TPCSector2DEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TPCSector2D*>(obj);

  fShowMax->SetState(fM->fShowMax ? kButtonDown : kButtonUp);
  SetupAverage();

  fUseTexture->SetState(fM->fUseTexture ? kButtonDown : kButtonUp);
}

/**************************************************************************/

void TPCSector2DEditor::DoShowMax()
{
  fM->SetShowMax(fShowMax->IsOn());
  SetupAverage();
  Update();
}

void TPCSector2DEditor::DoAverage()
{
  fM->SetAverage(fAverage->IsOn());
  Update();
}

void TPCSector2DEditor::SetupAverage()
{
  if(fM->fShowMax) {
    fAverage->SetEnabled(kFALSE);
  } else {
    fAverage->SetEnabled(kTRUE);
    fAverage->SetState(fM->fAverage ? kButtonDown : kButtonUp);
  }
}

/**************************************************************************/

void TPCSector2DEditor::DoUseTexture()
{
  fM->fUseTexture = fUseTexture->IsOn();
  Update();
}
