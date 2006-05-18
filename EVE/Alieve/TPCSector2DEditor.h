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
protected:
  TPCSector2D* fM; // fModel dynamic-casted to TPCSector2DEditor

  TGCheckButton*   fUseTexture;

  TGNumberEntry*   fSectorID;

  TGCheckButton*   fRnrInn;
  TGCheckButton*   fRnrOut1;
  TGCheckButton*   fRnrOut2;

  TGLabel*         fThresholdLabel;
  TGLabel*         fMaxValLabel;
  TGHSlider*       fThreshold;
  TGHSlider*       fMaxVal;

  TGCheckButton*   fShowMax;
  TGCheckButton*   fAverage;

  TGNumberEntry*   fMinTime;
  TGNumberEntry*   fMaxTime;
  TGDoubleHSlider* fTime;

public:
  TPCSector2DEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TPCSector2DEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);
  virtual Bool_t CanEditMainColor()  { return true; }

  // void DoXYZZ();
  void DoUseTexture();

  void DoSectorID();

  void DoRnrInn();
  void DoRnrOut1();
  void DoRnrOut2();

  void DoThreshold();
  void DoMaxVal();
  void DoShowMax();
  void DoAverage();
  void SetupAverage();
  void DoMinTime();
  void DoMaxTime();
  void DoTime();
 
  ClassDef(TPCSector2DEditor, 0); // Editor for TPCSector2D
}; // endclass TPCSector2DEditor

}

#endif
