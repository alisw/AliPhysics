// $Header$

#ifndef REVE_NLTProjectorEditor_H
#define REVE_NLTProjectorEditor_H

#include <TGedFrame.h>

class TGComboBox;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class NLTProjector;
class RGValuator;

class NLTProjectorEditor : public TGedFrame
{
private:
  NLTProjectorEditor(const NLTProjectorEditor&);            // Not implemented
  NLTProjectorEditor& operator=(const NLTProjectorEditor&); // Not implemented

protected:
  NLTProjector  *fM; // fModel dynamic-casted to NLTProjectorEditor

  TGComboBox    *fType;
  RGValuator    *fDistortion;
  RGValuator    *fFixedRadius;
  RGValuator    *fCurrentDepth;

  // axis
  TGColorSelect *fAxisColor;
  TGComboBox    *fSIMode; 
  TGNumberEntry *fSILevel;

public:
  NLTProjectorEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~NLTProjectorEditor();

  virtual void SetModel(TObject* obj);

  void DoSplitInfoMode(Int_t type);
  void DoSplitInfoLevel();
  void DoAxisColor(Pixel_t pixel);

  void DoType(Int_t type);
  void DoDistortion();
  void DoFixedRadius();  
  void DoCurrentDepth();

  ClassDef(NLTProjectorEditor, 0); // Editor for NLTProjector
}; // endclass NLTProjectorEditor

}

#endif
