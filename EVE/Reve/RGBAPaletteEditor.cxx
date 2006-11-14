// $Header$

#include "RGBAPaletteEditor.h"
#include <Reve/RGBAPalette.h>
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


RGBAPaletteSubEditor::RGBAPaletteSubEditor(const TGWindow* p) :
  TGVerticalFrame(p)
{
  // Create weeds.

  Int_t labelW = 42;

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

    fWrap = new TGCheckButton(f, "Wrap");
    f->AddFrame(fWrap, new TGLayoutHints(kLHintsTop, 3, 1, 1, 0));
    fWrap->Connect("Toggled(Bool_t)",
		   "Reve::RGBAPaletteSubEditor", this, "DoWrap()");
  }

  fMinVal = new RGValuator(this, "MinVal", 200, 0);
  fMinVal->SetNELength(4);
  fMinVal->SetLabelWidth(labelW);
  fMinVal->Build();
  fMinVal->GetSlider()->SetWidth(120);
  fMinVal->Connect("ValueSet(Double_t)",
		   "Reve::RGBAPaletteSubEditor", this, "DoMinVal()");
  AddFrame(fMinVal, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

  fMaxVal = new RGValuator(this,"MaxVal", 200, 0);
  fMaxVal->SetNELength(4);
  fMaxVal->SetLabelWidth(labelW);
  fMaxVal->Build();
  fMaxVal->GetSlider()->SetWidth(120);
  fMaxVal->SetLimits(0, 500);
  fMaxVal->Connect("ValueSet(Double_t)",
		   "Reve::RGBAPaletteSubEditor", this, "DoMaxVal()");
  AddFrame(fMaxVal, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));

}

/**************************************************************************/

void RGBAPaletteSubEditor::SetModel(RGBAPalette* p)
{
  fM = p;

  fMinVal->SetValue(fM->fMinVal);
  fMaxVal->SetValue(fM->fMaxVal);
  fMaxVal->SetLimits(fM->fLowLimit, fM->fHighLimit);
  fMinVal->SetLimits(fM->fLowLimit, fM->fHighLimit);
  fInterpolate->SetState(fM->fInterpolate ? kButtonDown : kButtonUp);
  fWrap->SetState(fM->fWrap ? kButtonDown : kButtonUp);
  fDefaultColor->SetColor(TColor::Number2Pixel(fM->GetDefaultColor()), kFALSE);
}

/**************************************************************************/

void RGBAPaletteSubEditor::Changed()
{
  Emit("Changed()");
}

/**************************************************************************/

void RGBAPaletteSubEditor::DoMinVal()
{
  fM->SetMin((Int_t) fMinVal->GetValue());
  fMinVal->SetValue(fM->fMinVal);
  Changed();
}

void RGBAPaletteSubEditor::DoMaxVal()
{
  fM->SetMax((Int_t) fMaxVal->GetValue());
  fMaxVal->SetValue(fM->fMaxVal);
  Changed();
}


void RGBAPaletteSubEditor::DoInterpolate()
{
  fM->SetInterpolate(fInterpolate->IsOn());
  Changed();
}

void RGBAPaletteSubEditor::DoWrap()
{
  fM->SetWrap(fWrap->IsOn());
  Changed();
}

void RGBAPaletteSubEditor::DoDefaultColor(Pixel_t color)
{
  fM->SetDefaultColor(color);
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
