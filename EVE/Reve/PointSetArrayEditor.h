// $Header$

#ifndef REVE_PointSetArrayEditor_H
#define REVE_PointSetArrayEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;

namespace Reve {

class PointSetArray;

class PointSetArrayEditor : public TGedFrame
{
protected:
  PointSetArray* fM; // fModel dynamic-casted to PointSetArrayEditor

  TGDoubleHSlider* fSlider;

public:
  PointSetArrayEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~PointSetArrayEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoScroll();

  ClassDef(PointSetArrayEditor, 1); // Editor for PointSetArray
}; // endclass PointSetArrayEditor

}

#endif
