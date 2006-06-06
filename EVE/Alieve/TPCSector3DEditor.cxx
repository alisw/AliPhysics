// $Header$

#include "TPCSector3DEditor.h"
#include <Alieve/TPCSector3D.h>

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
// TPCSector3DEditor
//

ClassImp(TPCSector3DEditor)

TPCSector3DEditor::TPCSector3DEditor(const TGWindow *p, Int_t id, Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("TPCSector3D");

  fRnrFrame = new TGCheckButton(this, "ShowFrame");
  AddFrame(fRnrFrame, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrFrame->Connect
    ("Toggled(Bool_t)","Alieve::TPCSector3DEditor", this, "DoRnrFrame()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Drift Velocity factor:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 25, 2, 1, 1));
    fDriftVel = new TGNumberEntry(f, 0., 6, -1,
                       TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAPositive,
                       TGNumberFormat::kNELLimitMinMax, 0.001, 1000.0);
    fDriftVel->GetNumberEntry()->SetToolTipText("Drift velocity factor.");
    f->AddFrame(fDriftVel, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fDriftVel->Associate(f);
    fDriftVel->Connect("ValueSet(Long_t)", "Alieve::TPCSector3DEditor", this, "DoDriftVel()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  // Register the editor.
  TClass *cl = TPCSector3D::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

TPCSector3DEditor::~TPCSector3DEditor()
{}

/**************************************************************************/

void TPCSector3DEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t /*event*/)
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(TPCSector3D::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<TPCSector3D*>(fModel);

  fRnrFrame->SetState(fM->fRnrFrame ? kButtonDown : kButtonUp);
  fDriftVel->SetNumber(fM->fDriftVel);

  SetActive();
}

/**************************************************************************/

void TPCSector3DEditor::DoRnrFrame()
{
  fM->SetRnrFrame(fRnrFrame->IsOn());
  Update();
}

void TPCSector3DEditor::DoDriftVel()
{
  fM->SetDriftVel((Float_t) fDriftVel->GetNumber());
  Update();
}
