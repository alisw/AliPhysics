// $Header$

#ifndef ALIEVE_TOFDigitsInfoEditor_H
#define ALIEVE_TOFDigitsInfoEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class TOFDigitsInfo;

class TOFDigitsInfoEditor : public TGedFrame
{
private:
  TOFDigitsInfoEditor(const TOFDigitsInfoEditor&);            // Not implemented
  TOFDigitsInfoEditor& operator=(const TOFDigitsInfoEditor&); // Not implemented

protected:
  TOFDigitsInfo* fM; // fModel dynamic-casted to TOFDigitsInfoEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;

public:
  TOFDigitsInfoEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TOFDigitsInfoEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(TOFDigitsInfoEditor, 0); // Editor for TOFDigitsInfo
}; // endclass TOFDigitsInfoEditor

}

#endif
