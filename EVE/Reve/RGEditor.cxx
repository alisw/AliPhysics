// $Header$

#include "RGEditor.h"
#include "RenderElement.h"
#include "ReveManager.h"

#include <TGedFrame.h>
#include <TGCanvas.h>
#include <TCanvas.h>

//______________________________________________________________________
// RGEditor
//

using namespace Reve;

ClassImp(RGEditor)

RGEditor::RGEditor(TCanvas* canvas, Int_t width, Int_t height) :
  TGedEditor(canvas),
  fRnrElement(0),
  fObject    (0)
{
  Resize(width, height);

  // Fix priority for TAttMarkerEditor.
  TClass* amClass = TClass::GetClass("TAttMarker");
  TClass* edClass = TClass::GetClass("TAttMarkerEditor");
  TGWindow *exroot = (TGWindow*) fClient->GetRoot();
  fClient->SetRoot(fTabContainer);
  SetFrameCreator(this);
  TGedFrame *frame = reinterpret_cast<TGedFrame*>(edClass->New());
  frame->SetModelClass(amClass);
  {
    Int_t off = edClass->GetDataMemberOffset("fPriority");
    if(off == 0)
      printf("ojej!\n");
    else
      * (Int_t*) (((char*)frame) + off) = 1;
  }
  SetFrameCreator(0);
  fClient->SetRoot(exroot);
  fFrameMap.Add(amClass, frame);
}

RenderElement* RGEditor::GetRnrElement() const
{
  return (fModel == fObject) ? fRnrElement : 0;
}

void RGEditor::DisplayRenderElement(RenderElement* re)
{
  fRnrElement = re;
  fObject     = fRnrElement ? fRnrElement->GetEditorObject() : 0;
  TGedEditor::SetModel(fPad, fObject, kButton1Down);
}

void RGEditor::DisplayObject(TObject* obj)
{
  fRnrElement = dynamic_cast<RenderElement*>(obj);
  fObject     = obj;
  TGedEditor::SetModel(fPad, obj, kButton1Down);
}

/**************************************************************************/

void RGEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t event)
{
  // !!!! do something so that such calls from elswhere will also
  // now the render element

  fRnrElement = dynamic_cast<RenderElement*>(obj);
  fObject     = obj;
  TGedEditor::SetModel(pad, obj, event);
}

void RGEditor::Update(TGedFrame* /*gframe*/)
{
  // Virtual method from TGedEditor ... called on every change.

  if (fRnrElement) {
    fRnrElement->UpdateItems();
    fRnrElement->ElementChanged();
  }

  gReve->Redraw3D();
}

/**************************************************************************/

/*
// Attempt to enable mouse-wheel in geditor -- failed.
Bool_t RGEditor::HandleButton(Event_t *event)
{
  // Handle mouse button event in container.

  printf("odfjgsf\n");
  if (event->fCode == kButton4 || event->fCode == kButton5) {
    return fCan->GetContainer()->HandleButton(event);
  } else {
    return TGedEditor::HandleButton(event);
  }
}
*/
