// $Header$

#ifndef REVE_GridStepperEditor_H
#define REVE_GridStepperEditor_H

#include <TGedFrame.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class GridStepper;
class RGValuator;

class GridStepperSubEditor : public TGVerticalFrame
{
private:
  GridStepperSubEditor(const GridStepperSubEditor&);            // Not implemented
  GridStepperSubEditor& operator=(const GridStepperSubEditor&); // Not implemented

protected:
  GridStepper   *fM;

  RGValuator    *fNx;
  RGValuator    *fNy;
  RGValuator    *fNz;    
  RGValuator    *fDx; 
  RGValuator    *fDy; 
  RGValuator    *fDz;

public:
  GridStepperSubEditor(const TGWindow* p);
  virtual ~GridStepperSubEditor() {}

  void SetModel(GridStepper* m);

  void Changed(); //*SIGNAL*

  void DoNs();
  void DoDs();

  ClassDef(GridStepperSubEditor, 0) // Sub-editor for GridStepper
};


class GridStepperEditor : public TGedFrame
{
private:
  GridStepperEditor(const GridStepperEditor&);            // Not implemented
  GridStepperEditor& operator=(const GridStepperEditor&); // Not implemented

protected:
  GridStepper            *fM;  // fModel dynamic-casted to GridStepper
  GridStepperSubEditor   *fSE;

public:
  GridStepperEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~GridStepperEditor() {}

  virtual void SetModel(TObject* obj);

  ClassDef(GridStepperEditor, 0) // Editor for GridStepper
};

}

#endif
