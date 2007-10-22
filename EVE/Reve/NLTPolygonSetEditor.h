#ifndef REVE_NLTPolygonSetEditor_H
#define REVE_NLTPolygonSetEditor_H

#include <TGedFrame.h>
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class NLTPolygonSet;
class NLTPolygonSetEditor : public TGedFrame
{
  NLTPolygonSetEditor(const NLTPolygonSetEditor&);            // Not implemented
  NLTPolygonSetEditor& operator=(const NLTPolygonSetEditor&); // Not implemented

protected:
  NLTPolygonSet* fPS; // fModel dynamic-casted to NLTPolygonSetEditor

  // TGColorSelect       *fFillColor;      // fill color widget
  
  TGNumberEntry       *fLineWidth;
  TGColorSelect       *fLineColor;      // fill color widget

public:
  NLTPolygonSetEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		   UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~NLTPolygonSetEditor();

  virtual void SetModel(TObject* obj);

  //virtual void DoFillColor(Pixel_t color);

  virtual void DoLineWidth();
  virtual void DoLineColor(Pixel_t color);

  ClassDef(NLTPolygonSetEditor, 0); // Editor for NLTPolygonSet
}; // endclass  NLTPolygonSetEditor
}
#endif
