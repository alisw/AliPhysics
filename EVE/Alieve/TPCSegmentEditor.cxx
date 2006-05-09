// $Header$

#include "TPCSegmentEditor.h"
#include <Alieve/TPCSegment.h>

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
// TPCSegmentEditor
//

ClassImp(TPCSegmentEditor)

TPCSegmentEditor::TPCSegmentEditor(const TGWindow *p, Int_t id, Int_t width, Int_t height,
				   UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("TPCSegment");

  //!!! create the widgets here ...

  // Register the editor.
  // TClass *cl = TPCSegment::Class();
  //  TGedElement *ge = new TGedElement;
  // ge->fGedFrame = this;
  //  ge->fCanvas = 0;
  //  cl->GetEditorList()->Add(ge);

  fUseTexture = new TGCheckButton(this, "UseTexture");
  AddFrame(fUseTexture, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fUseTexture->Connect
    ("Toggled(Bool_t)","Alieve::TPCSegmentEditor", this, "DoUseTexture()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "SegmentID:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 15, 2, 1, 1));

    fSegmentID = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
				   TGNumberFormat::kNELLimitMinMax, 0, 35);

    fSegmentID->GetNumberEntry()->SetToolTipText("0-18 front plate 18-36 back plate");
    f->AddFrame(fSegmentID, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fSegmentID->Associate(this);
    fSegmentID->Connect("ValueSet(Long_t)",
			"Alieve::TPCSegmentEditor", this, "DoSegmentID()");		    

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));  
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Treshold:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 15, 2, 1, 1));

    fTreshold = new TGHSlider(f, 150);
    fTreshold->SetRange(0,10);
    fTreshold->Associate(this);
    f->AddFrame(fTreshold, new TGLayoutHints(kLHintsLeft, 0, 5));
    fTreshold->Connect("PositionChanged(Int_t)",
		       "Alieve::TPCSegmentEditor", this, "DoTreshold()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "MaxValue:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 2, 1, 1));
    fMaxVal = new TGHSlider(f, 150);
    fMaxVal->SetRange(0,100);
    fMaxVal->Associate(this);
    f->AddFrame(fMaxVal, new TGLayoutHints(kLHintsLeft, 0, 5));
    fMaxVal->Connect("PositionChanged(Int_t)",
		     "Alieve::TPCSegmentEditor", this, "DoMaxVal()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  fShowMax = new TGCheckButton(this, "ShowMax");
  AddFrame(fShowMax, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
  fShowMax->Connect("Toggled(Bool_t)","Alieve::TPCSegmentEditor", this, "DoShowMax()");

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel *l = new TGLabel(f, "Time Range:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 1, 2, 1, 1));

    fTime = new TGDoubleHSlider(f);
    fTime->SetRange(0, 500);
    fTime->Resize(160, 20);
    f->AddFrame(fTime);//, new TGLayoutHints(kLHintsLeft, 0, 5));
    fTime->Connect("PositionChanged()", "Alieve::TPCSegmentEditor",
		   this, "DoTime()");
    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));
  }
  // What is this crap?
  TClass *cl = TPCSegmentEditor::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

TPCSegmentEditor::~TPCSegmentEditor()
{}

/**************************************************************************/

void TPCSegmentEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t )
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(TPCSegment::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<TPCSegment*>(fModel);

  fUseTexture->SetState(fM->fUseTexture ? kButtonDown : kButtonUp);
  fSegmentID->SetNumber(fM->fID);
  fTreshold->SetPosition(fM->fTreshold);//(fTreshold->GetMaxPosition()- fTreshold->GetMinPosition() ));

  fMaxVal->SetPosition(fM->fMaxVal);
  fTime->SetPosition(fM->fMinTime, fM->fMaxTime);

  fShowMax->SetState(fM->fShowMax ? kButtonDown : kButtonUp);
  //  fTime->SetPosition(fM->fMaxTime);
  SetActive();
}

/**************************************************************************/

void TPCSegmentEditor::DoUseTexture()
{
  fM->fUseTexture = fUseTexture->IsOn();
  Update();
}

void TPCSegmentEditor::DoSegmentID()
{
  fM->SetSegmentID((Int_t) fSegmentID->GetNumber());
  Update();
}

void TPCSegmentEditor::DoTreshold()
{
  fM->SetTreshold((Short_t) fTreshold->GetPosition());
  printf("DoTreshold %d \n",  fM->fTreshold);
  Update();
}

void TPCSegmentEditor::DoMaxVal()
{
  fM->SetMaxVal((Int_t) fMaxVal->GetPosition());
  printf("DoMaxVal %d \n",fM->fMaxVal);
  Update();
}

void TPCSegmentEditor::DoShowMax()
{
  fM->SetShowMax(fShowMax->IsOn());
  Update();
}

void TPCSegmentEditor::DoTime()
{ 
  Double_t min = fTime->GetMinPosition(), max = fTime->GetMaxPosition();
  printf("hslidor min=%f max=%f\n", min, max);
  fM->SetMinTime((Int_t) min);
  fM->SetMaxTime((Int_t) max);
  Update();
}
