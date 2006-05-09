// $Header$

#ifndef ALIEVE_TPCSegmentEditor_H
#define ALIEVE_TPCSegmentEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

namespace Alieve {

  class TPCSegment;

  class TPCSegmentEditor : public TGedFrame
  {
  protected:
    TPCSegment* fM; // fModel dynamic-casted to TPCSegmentEditor

    TGCheckButton*   fUseTexture;

    TGNumberEntry*   fSegmentID;

    TGHSlider*       fTreshold;
    TGHSlider*       fMaxVal;

    TGCheckButton*   fShowMax;

    TGDoubleHSlider* fTime;

  public:
    TPCSegmentEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
    ~TPCSegmentEditor();

    virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);
    virtual Bool_t CanEditMainColor()  { return true; }

    // void DoXYZZ();
    void DoUseTexture();

    void DoSegmentID();

    void DoTreshold();
    void DoMaxVal();
    void DoShowMax();
    void DoTime();
 
    ClassDef(TPCSegmentEditor, 1); // Editor for TPCSegment
  }; // endclass TPCSegmentEditor

}

#endif
