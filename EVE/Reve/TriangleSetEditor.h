// $Header$

#ifndef REVE_TriangleSetEditor_H
#define REVE_TriangleSetEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class ZTransSubEditor;
class TriangleSet;

class TriangleSetEditor : public TGedFrame
{
private:
  TriangleSetEditor(const TriangleSetEditor&);            // Not implemented
  TriangleSetEditor& operator=(const TriangleSetEditor&); // Not implemented

protected:
  TriangleSet* fM; // fModel dynamic-casted to TriangleSetEditor

  ZTransSubEditor *fHMTrans;

public:
  TriangleSetEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TriangleSetEditor();

  virtual void SetModel(TObject* obj);

  ClassDef(TriangleSetEditor, 1); // Editor for TriangleSet
}; // endclass TriangleSetEditor

}

#endif
