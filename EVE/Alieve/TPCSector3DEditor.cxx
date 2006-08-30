// $Header$

#include "TPCSector3DEditor.h"
#include <Alieve/TPCSector3D.h>

#include <Reve/RGValuators.h>

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

TPCSector3DEditor::TPCSector3DEditor(const TGWindow *p, Int_t id,
                                     Int_t width, Int_t height,
                                     UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back),
  fM(0),
  fRnrFrame(0), fDriftVel(0), fPointFrac(0), fPointSize(0)
{
  MakeTitle("TPCSector3D");

  Int_t labelW = 60;

  fRnrFrame = new TGCheckButton(this, "ShowFrame");
  AddFrame(fRnrFrame, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fRnrFrame->Connect
    ("Toggled(Bool_t)","Alieve::TPCSector3DEditor", this, "DoRnrFrame()");

  fDriftVel = new RGValuator(this, "Vdrift fac", 110, 0);
  fDriftVel->SetLabelWidth(labelW);
  fDriftVel->SetShowSlider(kFALSE);
  fDriftVel->SetNELength(6);
  fDriftVel->Build();
  fDriftVel->SetLimits(0.1, 10, 1, TGNumberFormat::kNESRealThree);
  fDriftVel->SetToolTip("Drift velocity factor");
  fDriftVel->Connect("ValueSet(Double_t)",
		     "Alieve::TPCSector3DEditor", this, "DoDriftVel()");
  AddFrame(fDriftVel, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fPointFrac = new RGValuator(this,"Point frac", 200, 0);
  fPointFrac->SetLabelWidth(labelW);
  fPointFrac->SetNELength(4);
  fPointFrac->Build();
  fPointFrac->GetSlider()->SetWidth(101 + 16);
  fPointFrac->SetLimits(0.0, 1.0, 101);
  fPointFrac->SetToolTip("Fraction of signal range displayed as points");
  fPointFrac->Connect("ValueSet(Double_t)",
		      "Alieve::TPCSector3DEditor", this, "DoPointFrac()");
  AddFrame(fPointFrac, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fPointSize = new RGValuator(this,"Point size", 200, 0);
  fPointSize->SetLabelWidth(labelW);
  fPointSize->SetShowSlider(kFALSE);
  fPointSize->SetNELength(4);
  fPointSize->Build();
  //fPointSize->GetSlider()->SetWidth(101 + 16);
  fPointSize->SetLimits(0.1, 32.0, 1, TGNumberFormat::kNESRealOne);
  fPointSize->SetToolTip("Size of displayed points");
  fPointSize->Connect("ValueSet(Double_t)",
		      "Alieve::TPCSector3DEditor", this, "DoPointSize()");
  AddFrame(fPointSize, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

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
  fDriftVel->SetValue(fM->fDriftVel);

  fPointFrac->SetValue(fM->fPointFrac);
  fPointSize->SetValue(fM->fPointSize);

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
  fM->SetDriftVel(fDriftVel->GetValue());
  Update();
}

void TPCSector3DEditor::DoPointFrac()
{
  fM->SetPointFrac(fPointFrac->GetValue());
  Update();
}

void TPCSector3DEditor::DoPointSize()
{
  fM->SetPointSize(fPointSize->GetValue());
  Update();
}

