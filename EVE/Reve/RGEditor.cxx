// $Header$

#include "RGEditor.h"

#include <TGedFrame.h>
#include <TGCanvas.h>
#include <TCanvas.h>

//______________________________________________________________________
// RGEditor
//

using namespace Reve;

ClassImp(RGEditor)

RGEditor::RGEditor(TCanvas* canvas) : TGedEditor(canvas)
{}

void RGEditor::DisplayObject(TObject* obj)
{
  fModel = obj;

  if(obj) {
    if(obj->IsA() != fClass && !obj->IsA()->InheritsFrom(fClass)) {
      fClass = obj->IsA();
      GetEditors();
    }
  } else {
    fCan->UnmapWindow();
    return;
  }

  TGFrameElement *el;
  TIter next(fStyle->GetList());
  while ((el = (TGFrameElement *) next())) {
    if ((el->fFrame)->InheritsFrom(TGedFrame::Class()))
      ((TGedFrame *)(el->fFrame))->SetModel(fPad, fModel, 0);
  }
  fCan->MapWindow();
}
