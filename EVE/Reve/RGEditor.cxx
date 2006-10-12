// $Header$

#include "RGEditor.h"
#include "RenderElement.h"
#include "RGTopFrame.h"

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

void RGEditor::DisplayRenderElement(RenderElement* re)
{
  fRnrElement = re;
  TObject* obj = fRnrElement ? fRnrElement->GetObject() : 0;
  SetModel(fPad, obj, kButton1Down);
}

void RGEditor::DisplayObject(TObject* obj)
{
  fRnrElement = 0;
  SetModel(fPad, obj, kButton1Down);
}

void RGEditor::Update(TGedFrame* /*gframe*/)
{
  // Virtual method from TGedEditor ... called on every change.

  if (fRnrElement) {
    fRnrElement->UpdateItems();
  }

  gReve->Redraw3D();
}
