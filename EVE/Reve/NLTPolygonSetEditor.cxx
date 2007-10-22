#include "NLTPolygonSetEditor.h"
#include <Reve/NLTPolygonSet.h>

#include <TColor.h>

#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGColorDialog.h>

using namespace Reve;


ClassImp(NLTPolygonSetEditor)

  NLTPolygonSetEditor::NLTPolygonSetEditor(const TGWindow *p,
					   Int_t width, Int_t height,
					   UInt_t options, Pixel_t back) :
    TGedFrame(p, width, height, options | kVerticalFrame, back),
    fPS(0),
    //  fFillColor(0),
    fLineWidth(0),
    fLineColor(0)
{
  MakeTitle("NLTPolygonSet");
  /*
    {
    TGCompositeFrame *f1 = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);
    TGLabel *l = new TGLabel(f1, "FillColor:");
    f1->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 25, 2, 1, 1));
    fFillColor = new TGColorSelect(f1, 0, -1);
    fFillColor->Connect("ColorSelected(Pixel_t)", "Reve::NLTPolygonSetEditor", this, "DoFillColor(Pixel_t)");
    f1->AddFrame(fFillColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
    AddFrame(f1, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
    }
  */

  {
    TGCompositeFrame *f = new TGCompositeFrame(this, 80, 20, kHorizontalFrame);

    TGLabel *l = new TGLabel(f, "LineColor:");
    f->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 25, 2, 1, 1));
    fLineColor = new TGColorSelect(f, 0, -1);
    fLineColor->Connect("ColorSelected(Pixel_t)", "Reve::NLTPolygonSetEditor", this, "DoLineColor(Pixel_t)");
    f->AddFrame(fLineColor, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));
 
    fLineWidth = new TGNumberEntry(f, 0., 6, -1, 
				   TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
				   TGNumberFormat::kNELLimitMinMax, 0.1, 2000.0);
    fLineWidth->GetNumberEntry()->SetToolTipText("Line witdth of outline.");
    fLineWidth->Connect("ValueSet(Long_t)", "Reve::NLTPolygonSetEditor", this, "DoLineWidth()");
    f->AddFrame(fLineWidth, new TGLayoutHints(kLHintsLeft, 1, 1, 1, 1));

    AddFrame(f, new TGLayoutHints(kLHintsTop, 1, 1, 0, 0));
  }
}


/**************************************************************************/
NLTPolygonSetEditor::~NLTPolygonSetEditor()
{}

/**************************************************************************/
/*
void NLTPolygonSetEditor::DoFillColor(Pixel_t pixel)
{
  printf("do fill color \n");
  fPS->SetFillColor(pixel);
  Update();
}
*/
/**************************************************************************/
void NLTPolygonSetEditor::DoLineWidth()
{
  // Double_t lw = fLineWidth->GetNumber();
  //fPS->SetLineWidth(lw);
  fPS->fLineWidth = fLineWidth->GetNumber();
  Update();
}

/**************************************************************************/
void NLTPolygonSetEditor::DoLineColor(Pixel_t pixel)
{
  fPS->SetLineColor(pixel);
  Update();
}

/**************************************************************************/

void NLTPolygonSetEditor::SetModel(TObject* obj)
{
  fPS = dynamic_cast<NLTPolygonSet*>(obj);
  //  fFillColor->SetColor(TColor::Number2Pixel(fPS->GetFillColor()), kFALSE);
  fLineWidth->SetNumber(fPS->fLineWidth);
  fLineColor->SetColor(TColor::Number2Pixel(fPS->GetLineColor()), kFALSE);
}
