#ifndef ALIEVE_MUONChamberEditor_H
#define ALIEVE_MUONChamberEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

namespace Reve {

class RGValuator;

}

namespace Alieve {

class MUONChamber;

class MUONChamberEditor : public TGedFrame
{

  MUONChamberEditor(const MUONChamberEditor&);            // Not implemented
  MUONChamberEditor& operator=(const MUONChamberEditor&); // Not implemented
  
 protected:

  MUONChamber* fM; // fModel dynamic-casted to MUONChamberEditor

  Reve::RGValuator *fThreshold;   // digit ADC min
  Reve::RGValuator *fMaxVal;      // digit ADC max
  Reve::RGValuator *fClusterSize; // cluster point size
  Reve::RGValuator *fHitSize;     // hit point size

 public:

  MUONChamberEditor(const TGWindow* p = 0,
		    Int_t width = 170, Int_t height = 30, 
		    UInt_t options = kChildFrame, 
		    Pixel_t back = GetDefaultFrameBackground());
  
  virtual ~MUONChamberEditor();

  virtual void SetModel(TObject* obj);

  void DoThreshold();
  void DoMaxVal();
  void DoClusterSize();
  void DoHitSize();

  ClassDef(MUONChamberEditor, 0); // Editor for MUONChamber

};

}

#endif
