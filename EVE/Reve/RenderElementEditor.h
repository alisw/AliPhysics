// $Header$

#ifndef REVE_RenderElementEditor_H
#define REVE_RenderElementEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class RenderElement;
class ZTransSubEditor;

class RenderElementEditor : public TGedFrame
{
  RenderElementEditor(const RenderElementEditor&);            // Not implemented
  RenderElementEditor& operator=(const RenderElementEditor&); // Not implemented

protected:
  RenderElement* fRE; // fModel dynamic-casted to RenderElementEditor

  TGHorizontalFrame* fHFrame;
  TGCheckButton*     fRnrSelf;
  TGCheckButton*     fRnrChildren;
  TGColorSelect*     fMainColor;
  TGNumberEntry*     fTransparency;
  ZTransSubEditor*   fHMTrans;

public:
  RenderElementEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		      UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~RenderElementEditor();

  virtual void SetModel(TObject* obj);

  void DoRnrSelf();
  void DoRnrChildren();
  void DoMainColor(Pixel_t color);
  void DoTransparency();

   ClassDef(RenderElementEditor, 1); // Editor for RenderElement
}; // endclass RenderElementEditor

}

#endif
