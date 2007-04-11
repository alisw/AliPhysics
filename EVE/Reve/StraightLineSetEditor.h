// $Header$

#ifndef REVE_StraightLineSetEditor_H
#define REVE_StraightLineSetEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class StraightLineSet;

class StraightLineSetEditor : public TGedFrame
{
private:
  StraightLineSetEditor(const StraightLineSetEditor&);            // Not implemented
  StraightLineSetEditor& operator=(const StraightLineSetEditor&); // Not implemented

protected:
  StraightLineSet* fM; // fModel dynamic-casted to StraightLineSetEditor

  // Declare widgets
  TGCheckButton*     fRnrMarkers;
  TGCheckButton*     fRnrLines;

public:
  StraightLineSetEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~StraightLineSetEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoRnrMarkers();
  void DoRnrLines();

  ClassDef(StraightLineSetEditor, 1); // Editor for StraightLineSet
}; // endclass StraightLineSetEditor

}

#endif
