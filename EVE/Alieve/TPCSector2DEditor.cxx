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

    fRnrInn = new TGCheckButton(f, "Inner");
    f->AddFrame(fRnrInn, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrInn->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoRnrInn()");

    fRnrOut1 = new TGCheckButton(f, "Outer 1");
    f->AddFrame(fRnrOut1, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrOut1->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoRnrOut1()");

    fRnrOut2 = new TGCheckButton(f, "Outer 2");
    f->AddFrame(fRnrOut2, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrOut2->Connect("Toggled(Bool_t)","Alieve::TPCSector2DEditor", this, "DoRnrOut2()");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    fThresholdLabel = new TGLabel(f, "Threshold [XXX]:");
    f->AddFrame(fThresholdLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 2, 1, 1));

    fThreshold = new TGHSlider(f, 150);
    fThreshold->SetRange(0,149);
    fThreshold->Associate(this);
    f->AddFrame(fThreshold, new TGLayoutHints(kLHintsLeft, 0, 5));
    fThreshold->Connect("PositionChanged(Int_t)",
			"Alieve::TPCSector2DEditor", this, "DoThreshold()");
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

  {
    TGVerticalFrame* vf = new TGVerticalFrame(this);
    
    TGHorizontalFrame* f = new TGHorizontalFrame(vf);

    TGLabel *l = new TGLabel(f, "Time Range:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 2, 1, 1));

    fMinTime = new TGNumberEntry(f, 0., 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax, 0, 1023);
    fMinTime->GetNumberEntry()->SetToolTipText("MinTime");
    f->AddFrame(fMinTime, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMinTime->Associate(this);
    fMinTime->Connect("ValueSet(Long_t)", "Alieve::TPCSector2DEditor", this, "DoMinTime()");

    fMaxTime = new TGNumberEntry(f, 0., 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax, 0, 1023);
    fMaxTime->GetNumberEntry()->SetToolTipText("MaxTime");
    f->AddFrame(fMaxTime, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fMaxTime->Associate(this);
    fMaxTime->Connect("ValueSet(Long_t)", "Alieve::TPCSector2DEditor", this, "DoMaxTime()");

    vf->AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

    fTime = new TGDoubleHSlider(vf);
    fTime->SetRange(0, 1023);
    fTime->Resize(250, 20);
    vf->AddFrame(fTime);//, new TGLayoutHints(kLHintsLeft, 0, 5));
    fTime->Connect("PositionChanged()", "Alieve::TPCSector2DEditor",
		   this, "DoTime()");

    AddFrame(vf, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
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
  fRnrInn->SetState(fM->fRnrInn   ? kButtonDown : kButtonUp);
  fRnrOut1->SetState(fM->fRnrOut1 ? kButtonDown : kButtonUp);
  fRnrOut2->SetState(fM->fRnrOut2 ? kButtonDown : kButtonUp);
  fThresholdLabel->SetText(Form("Threshold [%3d]:", fM->fThreshold));
  fThreshold->SetPosition(fM->fThreshold);

  fMaxValLabel->SetText(Form("MaxValue [%3d]:", fM->fMaxVal));
  fMaxVal->SetPosition(fM->fMaxVal);
  fMinTime->SetNumber(fM->fMinTime);
  fMaxTime->SetNumber(fM->fMaxTime);
  fTime->SetPosition(fM->fMinTime, fM->fMaxTime);

  fShowMax->SetState(fM->fShowMax ? kButtonDown : kButtonUp);
  SetupAverage();

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

void TPCSector2DEditor::DoRnrInn()
{
  fM->SetRnrInn(fRnrInn->IsOn());
  Update();
}

void TPCSector2DEditor::DoRnrOut1()
{
  fM->SetRnrOut1(fRnrOut1->IsOn());
  Update();
}

void TPCSector2DEditor::DoRnrOut2()
{
  fM->SetRnrOut2(fRnrOut2->IsOn());
  Update();
}


void TPCSector2DEditor::DoThreshold()
{
  fM->SetThreshold((Short_t) fThreshold->GetPosition());
  fThresholdLabel->SetText(Form("Threshold [%3d]:", fM->fThreshold));
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

void TPCSector2DEditor::DoMinTime()
{
  Int_t minTime = (Int_t) fMinTime->GetNumber();
  if(minTime > fM->fMaxTime) {
    minTime = fM->fMaxTime;
    fMinTime->SetNumber(minTime);
  }
  fM->SetMinTime(minTime);
  fTime->SetPosition(minTime, fM->fMaxTime);
  Update();
}

void TPCSector2DEditor::DoMaxTime()
{
  Int_t maxTime = (Int_t) fMaxTime->GetNumber();
  if(maxTime < fM->fMinTime) {
    maxTime = fM->fMinTime;
    fMaxTime->SetNumber(maxTime);
  }
  fM->SetMaxTime(maxTime);
  fTime->SetPosition(fM->fMinTime, maxTime);
  Update();
}

void TPCSector2DEditor::DoTime()
{ 
  Int_t min = (Int_t) TMath::Nint(fTime->GetMinPosition());
  Int_t max = (Int_t) TMath::Nint(fTime->GetMaxPosition());
  fM->SetMinTime(min);
  fM->SetMaxTime(max);
  fMinTime->SetNumber(min);
  fMaxTime->SetNumber(max);
  Update();
}
