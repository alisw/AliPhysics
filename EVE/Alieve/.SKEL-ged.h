// $Header$

#ifndef ALIEVE_CLASS_H
#define ALIEVE_CLASS_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class STEM;

class CLASS : public TGedFrame
{
protected:
  STEM* fM; // fModel dynamic-casted to CLASS

public:
  CLASS(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~CLASS();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  // void DoXYZZ();

  ClassDef(CLASS, 0); // Editor for STEM
}; // endclass CLASS

}

#endif
