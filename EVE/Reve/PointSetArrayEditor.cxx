// $Header$

#include "PointSetArrayEditor.h"
#include <Reve/PointSet.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// PointSetArrayEditor
//

ClassImp(PointSetArrayEditor)

PointSetArrayEditor::PointSetArrayEditor(const TGWindow *p, Int_t id, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, id, width, height, options | kVerticalFrame, back)
{
  fM = 0;
  MakeTitle("PointSetArray");

  fSlider = new TGDoubleHSlider(this);
  fSlider->Resize(260, 20);
  AddFrame(fSlider, new TGLayoutHints(kLHintsLeft, 0, 5));
  fSlider->Connect("PositionChanged()", "Reve::PointSetArrayEditor",
		   this, "DoScroll()");

  // Register the editor.
  TClass *cl = PointSetArray::Class();
  TGedElement *ge = new TGedElement;
  ge->fGedFrame = this;
  ge->fCanvas = 0;
  cl->GetEditorList()->Add(ge);
}

PointSetArrayEditor::~PointSetArrayEditor()
{}

/**************************************************************************/

void PointSetArrayEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t )
{
  fModel = 0;
  fPad   = 0;

  if (!obj || !obj->InheritsFrom(PointSetArray::Class()) || obj->InheritsFrom(TVirtualPad::Class())) {
    SetActive(kFALSE);
    return;
  }

  fModel = obj;
  fPad   = pad;

  fM = dynamic_cast<PointSetArray*>(fModel);

  printf("FullRange(%f, %f) Selected(%f,%f)\n",
	 fM->GetMin(), fM->GetMax(), fM->GetCurMin(), fM->GetCurMax());
  fSlider->SetRange(fM->GetMin(), fM->GetMax());
  fSlider->SetPosition(fM->GetCurMin(), fM->GetCurMax());

  SetActive();
}

/**************************************************************************/

void PointSetArrayEditor::DoScroll()
{
  Double_t min = fSlider->GetMinPosition(), max = fSlider->GetMaxPosition();
  printf("PointSet range: min=%f max=%f\n", min, max);
  fM->SetRange(min, max);
}
