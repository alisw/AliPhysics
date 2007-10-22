// $Header$

#ifndef REVE_RGEditor_H
#define REVE_RGEditor_H

#include <TGedEditor.h>

namespace Reve {

class RenderElement;

class RGEditor : public TGedEditor
{
  RGEditor(const RGEditor&);            // Not implemented
  RGEditor& operator=(const RGEditor&); // Not implemented

protected:
  RenderElement *fRnrElement; // Cached rnr-el pointer
  TObject       *fObject;     // Cached tobj pointer

public:
  RGEditor(TCanvas* canvas=0, Int_t width=250, Int_t height=400);
  virtual ~RGEditor() {}

  RenderElement* GetRnrElement() const;

  void DisplayRenderElement(RenderElement* re);
  void DisplayObject(TObject* obj);

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);
  virtual void Update(TGedFrame* gframe=0);

  // virtual Bool_t HandleButton(Event_t *event);

  ClassDef(RGEditor, 0);
}; // endclass RGEditor

}

#endif
