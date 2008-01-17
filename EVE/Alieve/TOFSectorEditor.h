// $Header$

#ifndef ALIEVE_TOFSectorEditor_H
#define ALIEVE_TOFSectorEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;

class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;
class TEveTransSubEditor;

namespace Alieve {
  
  class TOFSector;
  
  class TOFSectorEditor : public TGedFrame
  {
    //private:
    TOFSectorEditor(const TOFSectorEditor&);            // Not implemented
    TOFSectorEditor& operator=(const TOFSectorEditor&); // Not implemented

  protected:
    TOFSector* fM; // fModel dynamic-casted to TOFSectorEditor
    
    TEveGValuator* fSectorID;
    
    TGCheckButton*    fAutoTrans;

    TGCheckButton**    fPlate;

    TGCheckButton*    fPlate0;
    TGCheckButton*    fPlate1;
    TGCheckButton*    fPlate2;
    TGCheckButton*    fPlate3;
    TGCheckButton*    fPlate4;

    TEveGValuator* fThreshold;
    TEveGValuator* fMaxVal;   


    // Declare widgets
    // TGSomeWidget*   fXYZZ;
    
  public:
    TOFSectorEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
    virtual ~TOFSectorEditor();
    
    virtual void SetModel(TObject* obj);
    void DoSectorID();
    void DoAutoTrans(); 
    void DoPlate0();
    void DoPlate1();
    void DoPlate2();
    void DoPlate3();
    void DoPlate4();

    void DoPlate(Int_t nPlate);
    void DoThreshold();
    void DoMaxVal();


    // Declare callback/slot methods
    // void DoXYZZ();
    
    ClassDef(TOFSectorEditor, 0); // Editor for TOFSector
  }; // endclass TOFSectorEditor
  
}

#endif
