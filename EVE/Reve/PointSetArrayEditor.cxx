// $Header$

#include "PointSetArrayEditor.h"
#include <Reve/PointSet.h>
#include <Reve/RGValuators.h>

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

PointSetArrayEditor::PointSetArrayEditor(const TGWindow *p,
					 Int_t width, Int_t height,
					 UInt_t options, Pixel_t back) :
  TGedFrame(p,width, height, options | kVerticalFrame, back),
  fM(0),
  fRange(0)
{
  fM = 0;
  MakeTitle("PointSetArray");

  fRange = new RGDoubleValuator(this,"Range", 200, 0);
  fRange->SetNELength(6);
  //fRange->SetLabelWidth(labelW);
  fRange->Build();
  fRange->GetSlider()->SetWidth(224);
  fRange->Connect("ValueSet()",
		 "Reve::PointSetArrayEditor", this, "DoRange()");
  AddFrame(fRange, new TGLayoutHints(kLHintsTop, 1, 1, 2, 1));
}

PointSetArrayEditor::~PointSetArrayEditor()
{}

/**************************************************************************/

void PointSetArrayEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<PointSetArray*>(obj);

  // printf("FullRange(%f, %f) Selected(%f,%f)\n",
  //        fM->GetMin(), fM->GetMax(), fM->GetCurMin(), fM->GetCurMax());

  fRange->SetLimits(fM->fMin, fM->fMax, TGNumberFormat::kNESRealTwo);
  fRange->SetValues(fM->fCurMin, fM->fCurMax);
}

/**************************************************************************/

void PointSetArrayEditor::DoRange()
{
  fM->SetRange(fRange->GetMin(), fRange->GetMax());
  Update();
}
