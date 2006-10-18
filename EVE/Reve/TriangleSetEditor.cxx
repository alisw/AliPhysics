// $Header$

#include "TriangleSetEditor.h"
#include <Reve/TriangleSet.h>
#include <Reve/ZTransEditor.h>

#include <TVirtualPad.h>
#include <TColor.h>

#include <TGLabel.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGColorSelect.h>
#include <TGDoubleSlider.h>

using namespace Reve;

//______________________________________________________________________
// TriangleSetEditor
//

ClassImp(TriangleSetEditor)

TriangleSetEditor::TriangleSetEditor(const TGWindow *p, Int_t width, Int_t height,
	     UInt_t options, Pixel_t back) :
  TGedFrame(p, width, height, options | kVerticalFrame, back),
  fM(0),
  fHMTrans(0)
{
  MakeTitle("TriangleSet");

  fHMTrans = new ZTransSubEditor(this);
  fHMTrans->Connect("UseTrans()",     "Reve::TriangleSetEditor", this, "Update()");
  fHMTrans->Connect("TransChanged()", "Reve::TriangleSetEditor", this, "Update()");
  AddFrame(fHMTrans, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 0, 0, 0));
}

TriangleSetEditor::~TriangleSetEditor()
{
  delete fHMTrans;
}

/**************************************************************************/

void TriangleSetEditor::SetModel(TObject* obj)
{
  fM = dynamic_cast<TriangleSet*>(obj);

  fHMTrans->SetDataFromTrans(&fM->fHMTrans);
}
