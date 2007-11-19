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
  GridStepper   *fM;    // Model object.

  RGValuator    *fNx;   // Number of slots along x.
  RGValuator    *fNy;   // Number of slots along y.
  RGValuator    *fNz;   // Number of slots along z.
  RGValuator    *fDx;   // Step in the x direction.
  RGValuator    *fDy;   // Step in the y direction.
  RGValuator    *fDz;   // Step in the z direction.

public:
  GridStepperSubEditor(const TGWindow* p);
  virtual ~GridStepperSubEditor() {}

  void SetModel(GridStepper* m);

  void Changed(); //*SIGNAL*

  void DoNs();
  void DoDs();

  ClassDef(GridStepperSubEditor, 0) // Sub-editor for GridStepper class.
};


class GridStepperEditor : public TGedFrame
{
private:
  GridStepperEditor(const GridStepperEditor&);            // Not implemented
  GridStepperEditor& operator=(const GridStepperEditor&); // Not implemented

protected:
  GridStepper            *fM;   // Model object.
  GridStepperSubEditor   *fSE;  // Sub-editor containg GUI controls.

public:
  GridStepperEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~GridStepperEditor() {}

  virtual void SetModel(TObject* obj);

  ClassDef(GridStepperEditor, 0) // Editor for GridStepper class.
};

}

#endif
