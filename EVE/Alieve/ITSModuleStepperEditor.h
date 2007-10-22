// $Header$

#ifndef ALIEVE_ITSModuleStepperEditor_H
#define ALIEVE_ITSModuleStepperEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve
{
class GridStepperSubEditor;
}

namespace Alieve {

class ITSModuleStepper;

class ITSModuleStepperEditor : public TGedFrame
{
private:
  ITSModuleStepperEditor(const ITSModuleStepperEditor&);            // Not implemented
  ITSModuleStepperEditor& operator=(const ITSModuleStepperEditor&); // Not implemented

protected:
  ITSModuleStepper*     fM; // fModel dynamic-casted to ITSModuleStepperEditor

  Reve::GridStepperSubEditor* fStepper;
   
public:
  ITSModuleStepperEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~ITSModuleStepperEditor();

  virtual void SetModel(TObject* obj);

  void                  UpdateStore();
  ClassDef(ITSModuleStepperEditor, 0); // Editor for ITSModuleStepper
}; // endclass ITSModuleStepperEditor

}

#endif
