// $Header$

#ifndef ALIEVE_TPCSector2DEditor_H
#define ALIEVE_TPCSector2DEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;


namespace Alieve {

class TPCSector2D;

class TPCSector2DEditor : public TGedFrame
{
  TPCSector2DEditor(const TPCSector2DEditor&);            // Not implemented
  TPCSector2DEditor& operator=(const TPCSector2DEditor&); // Not implemented

protected:
  TPCSector2D* fM; // fModel dynamic-casted to TPCSector2DEditor

  TGCheckButton*   fShowMax;
  TGCheckButton*   fAverage;

  TGCheckButton*   fUseTexture;

public:
  TPCSector2DEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TPCSector2DEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoShowMax();
  void DoAverage();
  void SetupAverage();

  void DoUseTexture();

  ClassDef(TPCSector2DEditor, 0); // Editor for TPCSector2D
}; // endclass TPCSector2DEditor

}

#endif
