// $Header$

#include "TPCSector2DEditor.h"
#include <Alieve/TPCSector2D.h>

#include <TGButton.h>
#include <TGComboBox.h>
#include <TGLabel.h>

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
  fShowMax(0), fAverage(0), fUseTexture(0), fPickEmpty(0), fPickMode(0)
{
  MakeTitle("TPCSector2D");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fShowMax = new TGCheckButton(f, "ShowMax");
    f->AddFrame(fShowMax, new TGLayoutHints(kLHintsLeft, 3, 16, 1, 0));
    fShowMax->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoShowMax()");
    fAverage = new TGCheckButton(f, "Average");
    f->AddFrame(fAverage, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 0));
    fAverage->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoAverage()");
    AddFrame(f);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fUseTexture = new TGCheckButton(f, "UseTexture");
    f->AddFrame(fUseTexture, new TGLayoutHints(kLHintsTop, 3, 9, 1, 0));
    fUseTexture->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoUseTexture()");
    fPickEmpty = new TGCheckButton(f, "PickEmpty");
    f->AddFrame(fPickEmpty, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fPickEmpty->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoPickEmpty()");
    AddFrame(f);
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "PickMode");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 10, 1, 2));
    fPickMode = new TGComboBox(f);
    fPickMode->AddEntry("Print", 0);
    fPickMode->AddEntry("1D histo", 1);
    fPickMode->AddEntry("2D histo", 2);
    TGListBox* lb = fPickMode->GetListBox();
    lb->Resize(lb->GetWidth(), 3*16);
    fPickMode->Resize(80, 20);
    fPickMode->Connect("Selected(Int_t)", "Alieve::TPCSector2DEditor", this, "DoPickMode(Int_t)");
    f->AddFrame(fPickMode, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
    AddFrame(f);
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
  fPickEmpty->SetState(fM->fPickEmpty ? kButtonDown : kButtonUp);
  fPickMode->Select(fM->fPickMode, kFALSE);
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

void TPCSector2DEditor::DoPickEmpty()
{
  fM->fPickEmpty = fPickEmpty->IsOn();
  // Update();
}

void TPCSector2DEditor::DoPickMode(Int_t mode)
{
  fM->fPickMode = mode;
  // Update();
}
