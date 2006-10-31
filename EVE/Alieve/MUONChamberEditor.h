#ifndef ALIEVE_MUONChamberEditor_H
#define ALIEVE_MUONChamberEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

namespace Alieve {

class MUONChamber;

class MUONChamberEditor : public TGedFrame
{

  MUONChamberEditor(const MUONChamberEditor&);            // Not implemented
  MUONChamberEditor& operator=(const MUONChamberEditor&); // Not implemented
  
 protected:

  MUONChamber* fM; // fModel dynamic-casted to MUONChamberEditor

 public:

  MUONChamberEditor(const TGWindow* p,
		    Int_t width = 170, Int_t height = 30, 
		    UInt_t options = kChildFrame, 
		    Pixel_t back = GetDefaultFrameBackground());
  
  ~MUONChamberEditor();

  virtual void SetModel(TObject* obj);

  ClassDef(MUONChamberEditor, 0); // Editor for MUONChamber

};

}

#endif
