// $Header$

#include "TPCSectorVizEditor.h"
#include <Alieve/TPCSectorViz.h>

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
// TPCSectorVizEditor
//

ClassImp(TPCSectorVizEditor)

TPCSectorVizEditor::TPCSectorVizEditor(const TGWindow *p, Int_t id, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("TPCSectorViz");

  Int_t labelW = 60;

  fSectorID = new RGValuator(this, "SectorID", 110, 0);
  fSectorID->SetLabelWidth(labelW);
  fSectorID->SetShowSlider(kFALSE);
  fSectorID->SetNELength(4);
  fSectorID->Build();
  fSectorID->SetLimits(0, 35);
  fSectorID->SetToolTip("0-17 +z plate; 18-35 -z plate");
  fSectorID->Connect("ValueSet(Double_t)",
		     "Alieve::TPCSectorVizEditor", this, "DoSectorID()");
  AddFrame(fSectorID, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);

    fRnrInn = new TGCheckButton(f, "Inner");
    f->AddFrame(fRnrInn, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrInn->Connect("Toggled(Bool_t)","Alieve::TPCSectorVizEditor", this, "DoRnrInn()");

    fRnrOut1 = new TGCheckButton(f, "Outer 1");
    f->AddFrame(fRnrOut1, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrOut1->Connect("Toggled(Bool_t)","Alieve::TPCSectorVizEditor", this, "DoRnrOut1()");

    fRnrOut2 = new TGCheckButton(f, "Outer 2");
    f->AddFrame(fRnrOut2, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fRnrOut2->Connect("Toggled(Bool_t)","Alieve::TPCSectorVizEditor", this, "DoRnrOut2()");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }

  fThreshold = new RGValuator(this, "Threshold", 200, 0);
  fThreshold->SetNELength(4);
  fThreshold->SetLabelWidth(labelW);
  fThreshold->Build();
  fThreshold->GetSlider()->SetWidth(120);
  fThreshold->SetLimits(0,149);
  fThreshold->Connect("ValueSet(Double_t)",
		      "Alieve::TPCSectorVizEditor", this, "DoThreshold()");
  AddFrame(fThreshold, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fMaxVal = new RGValuator(this,"MaxVal", 200, 0);
  fMaxVal->SetNELength(4);
  fMaxVal->SetLabelWidth(labelW);
  fMaxVal->Build();
  fMaxVal->GetSlider()->SetWidth(120);
  fMaxVal->SetLimits(0, 299);
  fMaxVal->Connect("ValueSet(Double_t)",
		   "Alieve::TPCSectorVizEditor", this, "DoMaxVal()");
  AddFrame(fMaxVal, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fTime = new RGDoubleValuator(this,"Time", 200, 0);
  fTime->SetNELength(4);
  fTime->SetLabelWidth(labelW);
  fTime->Build();
  fTime->GetSlider()->SetWidth(224);
  fTime->SetLimits(0, 1023, TGNumberFormat::kNESInteger);
  fTime->Connect("ValueSet()",
		 "Alieve::TPCSectorVizEditor", this, "DoTime()");
  AddFrame(fTime, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));


  // Register the editor.
  TClass *cl = TPCSectorViz::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

TPCSectorVizEditor::~TPCSectorVizEditor()
{}

/**************************************************************************/

void TPCSectorVizEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t /*event*/)
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(TPCSectorViz::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<TPCSectorViz*>(fModel);

  fSectorID->SetValue(fM->fSectorID);

  fRnrInn ->SetState(fM->fRnrInn  ? kButtonDown : kButtonUp);
  fRnrOut1->SetState(fM->fRnrOut1 ? kButtonDown : kButtonUp);
  fRnrOut2->SetState(fM->fRnrOut2 ? kButtonDown : kButtonUp);

  fThreshold->SetValue(fM->fThreshold);
  fMaxVal->SetValue(fM->fMaxVal);

  fTime->SetValues(fM->fMinTime, fM->fMaxTime);

  SetActive();
}

/**************************************************************************/

void TPCSectorVizEditor::DoSectorID()
{
  fM->SetSectorID((Int_t) fSectorID->GetValue());
  Update();
}

/**************************************************************************/

void TPCSectorVizEditor::DoRnrInn()
{
  fM->SetRnrInn(fRnrInn->IsOn());
  Update();
}

void TPCSectorVizEditor::DoRnrOut1()
{
  fM->SetRnrOut1(fRnrOut1->IsOn());
  Update();
}

void TPCSectorVizEditor::DoRnrOut2()
{
  fM->SetRnrOut2(fRnrOut2->IsOn());
  Update();
}

/**************************************************************************/

void TPCSectorVizEditor::DoThreshold()
{
  fM->SetThreshold((Short_t) fThreshold->GetValue());
  Update();
}

void TPCSectorVizEditor::DoMaxVal()
{
  fM->SetMaxVal((Int_t) fMaxVal->GetValue());
  Update();
}

/**************************************************************************/

void TPCSectorVizEditor::DoTime()
{ 
  fM->SetMinTime((Int_t) fTime->GetMin());
  fM->SetMaxTime((Int_t) fTime->GetMax());
  Update();
}
