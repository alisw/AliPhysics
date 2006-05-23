// $Header$

#ifndef ALIEVE_TPCSectorVizEditor_H
#define ALIEVE_TPCSectorVizEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;


namespace Alieve {

class TPCSectorViz;

class TPCSectorVizEditor : public TGedFrame
{
protected:
  TPCSectorViz* fM; // fModel dynamic-casted to TPCSectorVizEditor

  TGNumberEntry*   fSectorID;

  TGCheckButton*   fRnrInn;
  TGCheckButton*   fRnrOut1;
  TGCheckButton*   fRnrOut2;

  TGLabel*         fThresholdLabel;
  TGLabel*         fMaxValLabel;
  TGHSlider*       fThreshold;
  TGHSlider*       fMaxVal;

  TGNumberEntry*   fMinTime;
  TGNumberEntry*   fMaxTime;
  TGDoubleHSlider* fTime;

public:
  TPCSectorVizEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TPCSectorVizEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoSectorID();

  void DoRnrInn();
  void DoRnrOut1();
  void DoRnrOut2();

  void DoThreshold();
  void DoMaxVal();

  void DoMinTime();
  void DoMaxTime();
  void DoTime();
 
  ClassDef(TPCSectorVizEditor, 0); // Editor for TPCSectorViz
}; // endclass TPCSectorVizEditor

}

#endif
