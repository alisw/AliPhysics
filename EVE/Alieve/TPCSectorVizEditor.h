// $Header$

#ifndef ALIEVE_TPCSectorVizEditor_H
#define ALIEVE_TPCSectorVizEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

namespace Reve {
class RGValuator;
class RGDoubleValuator;
}

namespace Alieve {

class TPCSectorViz;

class TPCSectorVizEditor : public TGedFrame
{
  TPCSectorVizEditor(const TPCSectorVizEditor&);            // Not implemented
  TPCSectorVizEditor& operator=(const TPCSectorVizEditor&); // Not implemented

protected:
  TPCSectorViz* fM; // fModel dynamic-casted to TPCSectorVizEditor

  Reve::RGValuator* fSectorID;
  TGCheckButton*    fTrans;

  TGCheckButton*    fRnrInn;
  TGCheckButton*    fRnrOut1;
  TGCheckButton*    fRnrOut2;

  Reve::RGValuator* fThreshold;
  Reve::RGValuator* fMaxVal;   

  Reve::RGDoubleValuator* fTime;

public:
  TPCSectorVizEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TPCSectorVizEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoSectorID();
  void DoTrans();

  void DoRnrInn();
  void DoRnrOut1();
  void DoRnrOut2();

  void DoThreshold();
  void DoMaxVal();

  void DoTime();
 
  ClassDef(TPCSectorVizEditor, 0); // Editor for TPCSectorViz
}; // endclass TPCSectorVizEditor

}

#endif
