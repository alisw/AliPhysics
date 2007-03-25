// $Header$

#ifndef ALIEVE_ITSModuleStepperGL_H
#define ALIEVE_ITSModuleStepperGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Alieve {

class ITSModuleStepper;

class ITSModuleStepperGL : public TGLObject
{
private:
  ITSModuleStepperGL(const ITSModuleStepperGL&);            // Not implemented
  ITSModuleStepperGL& operator=(const ITSModuleStepperGL&); // Not implemented

  void   RenderTriangle(Float_t dx, Float_t dy, Int_t id) const;
protected:
  ITSModuleStepper* fM; // fModel dynamic-casted to ITSModuleStepperGL

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  ITSModuleStepperGL();
  virtual ~ITSModuleStepperGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  // To support two-level selection
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(ITSModuleStepperGL, 0);
}; // endclass ITSModuleStepperGL

}

#endif
