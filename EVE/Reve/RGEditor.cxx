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
  SetModel(fPad, obj, kButton1Down);
}

void RGEditor::Update(TGedFrame* /*gframe*/)
{
  // Copy from TGeEditor ... need to do something better now.
  if (fPad) {
      fPad->Modified();
      fPad->Update();
   }

}
