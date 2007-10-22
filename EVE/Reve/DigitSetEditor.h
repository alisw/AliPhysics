// $Header$

#ifndef REVE_DigitSetEditor_H
#define REVE_DigitSetEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class DigitSet;

class RGValuator;
class RGDoubleValuator;
class ZTransSubEditor;

// It would be better to have button to change model to the palette
// object itself.
class RGBAPaletteSubEditor;

class DigitSetEditor : public TGedFrame
{
private:
  DigitSetEditor(const DigitSetEditor&);            // Not implemented
  DigitSetEditor& operator=(const DigitSetEditor&); // Not implemented

  void CreateInfoTab();
protected:
  DigitSet*              fM; // fModel dynamic-casted to DigitSetEditor

  ZTransSubEditor*      fHMTrans;
  RGBAPaletteSubEditor* fPalette;

  TGHorizontalFrame*    fHistoButtFrame;
  TGVerticalFrame*      fInfoFrame;

public:
  DigitSetEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~DigitSetEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoHisto();
  void DoRangeHisto();
  void PlotHisto(Int_t min, Int_t max);

  ClassDef(DigitSetEditor, 1); // Editor for DigitSet
}; // endclass DigitSetEditor

}

#endif
