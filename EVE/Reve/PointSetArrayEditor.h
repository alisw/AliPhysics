// $Header$

#ifndef REVE_PointSetArrayEditor_H
#define REVE_PointSetArrayEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class RGValuator;
class RGDoubleValuator;

class PointSetArray;

class PointSetArrayEditor : public TGedFrame
{
  PointSetArrayEditor(const PointSetArrayEditor&);            // Not implemented
  PointSetArrayEditor& operator=(const PointSetArrayEditor&); // Not implemented

protected:
  PointSetArray* fM; // fModel dynamic-casted to PointSetArrayEditor

  Reve::RGDoubleValuator* fRange;

public:
  PointSetArrayEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		      UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~PointSetArrayEditor();

  virtual void SetModel(TObject* obj);

  void DoRange();

  ClassDef(PointSetArrayEditor, 1); // Editor for PointSetArray
}; // endclass PointSetArrayEditor

}

#endif
