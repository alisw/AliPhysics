// $Header$

#ifndef REVE_CLASS_H
#define REVE_CLASS_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class STEM;

class CLASS : public TGedFrame
{
private:
  CLASS(const CLASS&);            // Not implemented
  CLASS& operator=(const CLASS&); // Not implemented

protected:
  STEM* fM; // fModel dynamic-casted to CLASS

  // Declare widgets
  // TGSomeWidget*   fXYZZ;

public:
  CLASS(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~CLASS();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(CLASS, 1); // Editor for STEM
}; // endclass CLASS

}

#endif
