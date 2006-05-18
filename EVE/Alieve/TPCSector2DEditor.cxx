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

TPCSector2DEditor::TPCSector2DEditor(const TGWindow *p, Int_t id, Int_t width, Int_t height,
				   UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("TPCSector2D");

  //!!! create the widgets here ...

  // Register the editor.
  // TClass *cl = TPCSector2D::Class();
  //  TGedElement *ge = new TGedElement;
  // ge->fGedFrame = this;
  //  ge->fCanvas = 0;
  //  cl->GetEditorList()->Add(ge);

  fUseTexture = new TGCheckButton(this, "UseTexture");
  AddFrame(fUseTexture, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fUseTexture->Connect
    ("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoUseTexture()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "SectorID:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 2, 1, 1));

    fSectorID = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				   TGNumberFormat::kNELLimitMinMax, 0, 35);

    fSectorID->GetNumberEntry()->SetToolTipText("0-17 +z plate 18-35 -z plate");
    f->AddFrame(fSectorID, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fSectorID->Associate(this);
    fSectorID->Connect("ValueSet(Long_t)",
			"Alieve::TPCSector2DEditor", this, "DoSectorID()");		    

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));  
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fThresholdLabel = new TGLabel(f, "threshold [XXX]:");
    f->AddFrame(fThresholdLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 2, 1, 1));

    fthreshold = new TGHSlider(f, 150);
    fthreshold->SetRange(0,149);
    fthreshold->Associate(this);
    f->AddFrame(fthreshold, new TGLayoutHints(kLHintsLeft, 0, 5));
    fthreshold->Connect("PositionChanged(Int_t)",
		       "Alieve::TPCSector2DEditor", this, "Dothreshold()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fMaxValLabel = new TGLabel(f, "MaxValue [XXX]:");
    f->AddFrame(fMaxValLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 2, 1, 1));
    fMaxVal = new TGHSlider(f, 150);
    fMaxVal->SetRange(0,299);
    fMaxVal->Associate(this);
    f->AddFrame(fMaxVal, new TGLayoutHints(kLHintsLeft, 0, 5));
    fMaxVal->Connect("PositionChanged(Int_t)",
		     "Alieve::TPCSector2DEditor", this, "DoMaxVal()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  fShowMax = new TGCheckButton(this, "ShowMax");
  AddFrame(fShowMax, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fShowMax->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoShowMax()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Time Range:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 2, 1, 1));

    fTime = new TGDoubleHSlider(f);
    fTime->SetRange(0, 500);
    fTime->Resize(160, 20);
    f->AddFrame(fTime);//, new TGLayoutHints(kLHintsLeft, 0, 5));
    fTime->Connect("PositionChanged()", "Alieve::TPCSector2DEditor",
		   this, "DoTime()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  // What is this crap?
  TClass *cl = TPCSector2DEditor::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

TPCSector2DEditor::~TPCSector2DEditor()
{}

/**************************************************************************/

void TPCSector2DEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t )
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(TPCSector2D::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<TPCSector2D*>(fModel);

  fUseTexture->SetState(fM->fUseTexture ? kButtonDown : kButtonUp);
  fSectorID->SetNumber(fM->fSectorID);
  fThresholdLabel->SetText(Form("threshold [%3d]:", fM->fthreshold));
  fthreshold->SetPosition(fM->fthreshold);

  fMaxValLabel->SetText(Form("MaxValue [%3d]:", fM->fMaxVal));
  fMaxVal->SetPosition(fM->fMaxVal);
  fTime->SetPosition(fM->fMinTime, fM->fMaxTime);

  fShowMax->SetState(fM->fShowMax ? kButtonDown : kButtonUp);

  SetActive();
}

/**************************************************************************/

void TPCSector2DEditor::DoUseTexture()
{
  fM->fUseTexture = fUseTexture->IsOn();
  Update();
}

void TPCSector2DEditor::DoSectorID()
{
  fM->SetSectorID((Int_t) fSectorID->GetNumber());
  Update();
}

void TPCSector2DEditor::Dothreshold()
{
  fM->Setthreshold((Short_t) fthreshold->GetPosition());
  fThresholdLabel->SetText(Form("threshold [%3d]:", fM->fthreshold));
  Update();
}

void TPCSector2DEditor::DoMaxVal()
{
  fM->SetMaxVal((Int_t) fMaxVal->GetPosition());
  fMaxValLabel->SetText(Form("MaxValue [%3d]:", fM->fMaxVal));
  Update();
}

void TPCSector2DEditor::DoShowMax()
{
  fM->SetShowMax(fShowMax->IsOn());
  Update();
}

void TPCSector2DEditor::DoTime()
{ 
  Double_t min = fTime->GetMinPosition(), max = fTime->GetMaxPosition();
  printf("hslidor min=%f max=%f\n", min, max);
  fM->SetMinTime((Int_t) min);
  fM->SetMaxTime((Int_t) max);
  Update();
}
