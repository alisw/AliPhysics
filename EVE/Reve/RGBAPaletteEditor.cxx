// $Header$

#include "RGBAPaletteEditor.h"
#include <Reve/RGBAPalette.h>
#include <Reve/RGValuators.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGColorSelect.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

using namespace Reve;


RGBAPaletteSubEditor::RGBAPaletteSubEditor(const TGWindow* p) :
  TGVerticalFrame(p),

  fM(0),

  fUndershootAction (0),
  fUnderColor       (0), 
  fOvershootAction  (0),
  fOverColor        (0),

  fMinMax(0),

  fInterpolate(0),
  fShowDefValue(0), 
  fDefaultColor(0)
{
  // Int_t labelW = 42;

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);

    fDefaultColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fDefaultColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fDefaultColor->Connect("ColorSelected(Pixel_t)",
			   "Reve::RGBAPaletteSubEditor", this, "DoDefaultColor(Pixel_t)");

    fInterpolate = new TGCheckButton(f, "Interpolate");
    f->AddFrame(fInterpolate, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fInterpolate->Connect("Toggled(Bool_t)",
			  "Reve::RGBAPaletteSubEditor", this, "DoInterpolate()");

    fShowDefValue = new TGCheckButton(f, "ShowDefValue");
    f->AddFrame(fShowDefValue, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fShowDefValue->Connect("Toggled(Bool_t)",
		   "Reve::RGBAPaletteSubEditor", this, "DoShowDefValue()");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
  }

  { // Undershoot
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "Under-mode");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 10, 1, 2));
    fUndershootAction = new TGComboBox(f);
    fUndershootAction->AddEntry("Cut", 0);
    fUndershootAction->AddEntry("Mark", 1);
    fUndershootAction->AddEntry("Clip", 2);
    fUndershootAction->AddEntry("Wrap", 3);
    TGListBox* lb = fUndershootAction->GetListBox();
    lb->Resize(lb->GetWidth(), 4*16);
    fUndershootAction->Resize(80, 20);
    fUndershootAction->Connect("Selected(Int_t)", "Reve::RGBAPaletteSubEditor", this,
		       "DoUndershootAction(Int_t)");
    f->AddFrame(fUndershootAction, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

    fUnderColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fUnderColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fUnderColor->Connect("ColorSelected(Pixel_t)",
			   "Reve::RGBAPaletteSubEditor", this, "DoUnderColor(Pixel_t)");
  
    AddFrame(f);
  }

  { // Overshoot
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "Over-mode");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 10, 1, 2));
    fOvershootAction = new TGComboBox(f);
    fOvershootAction->AddEntry("Cut", 0);
    fOvershootAction->AddEntry("Mark", 1);
    fOvershootAction->AddEntry("Clip", 2);
    fOvershootAction->AddEntry("Wrap", 3);
    TGListBox* lb = fOvershootAction->GetListBox();
    lb->Resize(lb->GetWidth(), 4*16);
    fOvershootAction->Resize(80, 20);
    fOvershootAction->Connect("Selected(Int_t)", "Reve::RGBAPaletteSubEditor", this,
		       "DoOvershootAction(Int_t)");
    f->AddFrame(fOvershootAction, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

    fOverColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fOverColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    fOverColor->Connect("ColorSelected(Pixel_t)",
			   "Reve::RGBAPaletteSubEditor", this, "DoOverColor(Pixel_t)");
  
    AddFrame(f);
  }

  fMinMax = new RGDoubleValuator(this,"Min / Max", 200, 0);
  fMinMax->SetNELength(4);
  fMinMax->SetLabelWidth(64);
  fMinMax->Build();
  fMinMax->GetSlider()->SetWidth(224);
  fMinMax->SetLimits(0, 1023, TGNumberFormat::kNESInteger);
  fMinMax->Connect("ValueSet()",
		   "Reve::RGBAPaletteSubEditor", this, "DoMinMax()");
  AddFrame(fMinMax, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

}

/**************************************************************************/

void RGBAPaletteSubEditor::SetModel(RGBAPalette* p)
{
  fM = p;

  fMinMax->SetValues(fM->fMinVal, fM->fMaxVal);
  fMinMax->SetLimits(fM->fLowLimit, fM->fHighLimit);

  fInterpolate->SetState(fM->fInterpolate ? kButtonDown : kButtonUp);
  fShowDefValue->SetState(fM->fShowDefValue ? kButtonDown : kButtonUp);
  fDefaultColor->SetColor(TColor::Number2Pixel(fM->GetDefaultColor()), kFALSE);

  fUnderColor->SetColor(TColor::Number2Pixel(fM->GetUnderColor()), kFALSE);
  fOverColor->SetColor(TColor::Number2Pixel(fM->GetOverColor()), kFALSE);

  fUndershootAction->Select(fM->fUndershootAction, kFALSE);
  fOvershootAction->Select(fM->fOvershootAction, kFALSE);
}

/**************************************************************************/

void RGBAPaletteSubEditor::Changed()
{
  Emit("Changed()");
}

/**************************************************************************/

void RGBAPaletteSubEditor::DoMinMax()
{ 
  fM->SetMinMax((Int_t) fMinMax->GetMin(), (Int_t) fMinMax->GetMax());
  Changed();
}

/**************************************************************************/

void RGBAPaletteSubEditor::DoInterpolate()
{
  fM->SetInterpolate(fInterpolate->IsOn());
  Changed();
}

void RGBAPaletteSubEditor::DoShowDefValue()
{
  fM->SetShowDefValue(fShowDefValue->IsOn());
  Changed();
}

void RGBAPaletteSubEditor::DoDefaultColor(Pixel_t color)
{
  fM->SetDefaultColor(color);
  Changed();
}

void RGBAPaletteSubEditor::DoUnderColor(Pixel_t color)
{
  fM->SetUnderColor(color);
  Changed();
}

void RGBAPaletteSubEditor::DoOverColor(Pixel_t color)
{
  fM->SetOverColor(color);
  Changed();
}

void RGBAPaletteSubEditor::DoUndershootAction(Int_t mode)
{
  fM->SetUndershootAction(mode);
  Changed();
}

void RGBAPaletteSubEditor::DoOvershootAction(Int_t mode)
{
  fM->SetOvershootAction(mode);
  Changed();
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// RGBAPaletteEditor
//

ClassImp(RGBAPaletteEditor)

RGBAPaletteEditor::RGBAPaletteEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM (0),
  fSE(0)
{
  MakeTitle("RGBAPalette");

  fSE = new RGBAPaletteSubEditor(this);
  AddFrame(fSE, new TGLayoutHints(kLHintsTop, 2, 0, 2, 2));
  fSE->Connect("Changed()", "Reve::RGBAPaletteEditor", this, "Update()");
}

RGBAPaletteEditor::~RGBAPaletteEditor()
{}

/**************************************************************************/

void RGBAPaletteEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<RGBAPalette*>(obj);
  fSE->SetModel(fM);
}
