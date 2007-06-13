// $Header$

#ifndef ALIEVE_TOFStripEditor_H
#define ALIEVE_TOFStripEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class TOFStrip;

class TOFStripEditor : public TGedFrame
{
private:
  TOFStripEditor(const TOFStripEditor&);            // Not implemented
  TOFStripEditor& operator=(const TOFStripEditor&); // Not implemented

protected:
  TOFStrip* fM; // fModel dynamic-casted to TOFStripEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;

public:
  TOFStripEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TOFStripEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(TOFStripEditor, 0); // Editor for TOFStrip
}; // endclass TOFStripEditor

}

#endif
