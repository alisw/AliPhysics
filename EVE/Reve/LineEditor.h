// $Header$

#ifndef REVE_LineEditor_H
#define REVE_LineEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class Line;

class LineEditor : public TGedFrame
{
private:
  LineEditor(const LineEditor&);            // Not implemented
  LineEditor& operator=(const LineEditor&); // Not implemented

protected:
  Line* fM; // fModel dynamic-casted to LineEditor

  TGCheckButton     *fRnrLine;
  TGCheckButton     *fRnrPoints;

public:
  LineEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~LineEditor();

  virtual void SetModel(TObject* obj);

  void DoRnrLine();
  void DoRnrPoints();

  ClassDef(LineEditor, 1); // Editor for Line
}; // endclass LineEditor

}

#endif
