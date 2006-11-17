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

  fUnderflowAction (0),
  fUnderColor      (0), 
  fOverflowAction  (0),
  fOverColor       (0),

  fMinMax(0),

  fInterpolate(0),
  fShowDefValue(0), 
  fDefaultColor(0)
{
  // Int_t labelW = 42;

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);

    fInterpolate = new TGCheckButton(f, "Interpolate");
    f->AddFrame(fInterpolate, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 0));
    fInterpolate->Connect("Toggled(Bool_t)",
			  "Reve::RGBAPaletteSubEditor", this, "DoInterpolate()");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));
  }

  {
    TGHorizontalFrame* f = new TGHorizontalFrame(this);

    fShowDefValue = new TGCheckButton(f, "Show default value");
    f->AddFrame(fShowDefValue, new TGLayoutHints(kLHintsLeft, 3, 1, 1, 0));
    fShowDefValue->Connect("Toggled(Bool_t)",
		   "Reve::RGBAPaletteSubEditor", this, "DoShowDefValue()");

    fDefaultColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fDefaultColor, new TGLayoutHints(kLHintsLeft|kLHintsTop, 3, 1, 0, 2));
    fDefaultColor->Connect("ColorSelected(Pixel_t)",
			   "Reve::RGBAPaletteSubEditor", this, "DoDefaultColor(Pixel_t)");

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 1, 0));
  }

  { // Underflow
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "Underflow:");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 15, 1, 2));
    fUnderflowAction = new TGComboBox(f);
    fUnderflowAction->AddEntry("Cut", 0);
    fUnderflowAction->AddEntry("Mark", 1);
    fUnderflowAction->AddEntry("Clip", 2);
    fUnderflowAction->AddEntry("Wrap", 3);
    TGListBox* lb = fUnderflowAction->GetListBox();
    lb->Resize(lb->GetWidth(), 4*16);
    fUnderflowAction->Resize(59, 20);
    fUnderflowAction->Connect("Selected(Int_t)", "Reve::RGBAPaletteSubEditor", this,
			      "DoUnderflowAction(Int_t)");
    f->AddFrame(fUnderflowAction, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));

    fUnderColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fUnderColor, new TGLayoutHints(kLHintsLeft|kLHintsTop, 1, 1, 0, 2));
    fUnderColor->Connect("ColorSelected(Pixel_t)",
			 "Reve::RGBAPaletteSubEditor", this, "DoUnderColor(Pixel_t)");
  
    AddFrame(f);
  }

  { // Overflow
    TGHorizontalFrame* f = new TGHorizontalFrame(this);
    TGLabel* lab = new TGLabel(f, "Overflow:");
    f->AddFrame(lab, new TGLayoutHints(kLHintsLeft|kLHintsBottom, 1, 20, 1, 2));
    fOverflowAction = new TGComboBox(f);
    fOverflowAction->AddEntry("Cut", 0);
    fOverflowAction->AddEntry("Mark", 1);
    fOverflowAction->AddEntry("Clip", 2);
    fOverflowAction->AddEntry("Wrap", 3);
    TGListBox* lb = fOverflowAction->GetListBox();
    lb->Resize(lb->GetWidth(), 4*16);
    fOverflowAction->Resize(59, 20);
    fOverflowAction->Connect("Selected(Int_t)", "Reve::RGBAPaletteSubEditor", this,
			     "DoOverflowAction(Int_t)");
    f->AddFrame(fOverflowAction, new TGLayoutHints(kLHintsLeft, 1, 2, 1, 1));

    fOverColor = new TGColorSelect(f, 0, -1);
    f->AddFrame(fOverColor, new TGLayoutHints(kLHintsLeft|kLHintsTop, 1, 1, 0, 2));
    fOverColor->Connect("ColorSelected(Pixel_t)",
			"Reve::RGBAPaletteSubEditor", this, "DoOverColor(Pixel_t)");
  
    AddFrame(f);
  }

  fMinMax = new RGDoubleValuator(this,"Main range:", 200, 0);
  fMinMax->SetNELength(5);
  fMinMax->SetLabelWidth(74);
  fMinMax->Build();
  fMinMax->GetSlider()->SetWidth(224);
  fMinMax->SetLimits(0, 1023, TGNumberFormat::kNESInteger);
  fMinMax->Connect("ValueSet()",
		   "Reve::RGBAPaletteSubEditor", this, "DoMinMax()");
  AddFrame(fMinMax, new TGLayoutHints(kLHintsTop, 1, 1, 1, 1));

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

  fUnderflowAction->Select(fM->fUnderflowAction, kFALSE);
  fOverflowAction->Select(fM->fOverflowAction, kFALSE);
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

void RGBAPaletteSubEditor::DoUnderflowAction(Int_t mode)
{
  fM->SetUnderflowAction(mode);
  Changed();
}

void RGBAPaletteSubEditor::DoOverflowAction(Int_t mode)
{
  fM->SetOverflowAction(mode);
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
