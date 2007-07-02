// $Header$

#ifndef ALIEVE_ITSModuleStepperGL_H
#define ALIEVE_ITSModuleStepperGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;
class TString;

namespace Alieve {

class ITSModuleStepper;

class ITSModuleStepperGL : public TGLObject
{
private:
  ITSModuleStepperGL(const ITSModuleStepperGL&);            // Not implemented
  ITSModuleStepperGL& operator=(const ITSModuleStepperGL&); // Not implemented

  void   RenderSymbol(Float_t dx, Float_t dy, Int_t id) const;
  void   RenderString(TString text , Float_t dy, Bool_t trans = kTRUE) const;
  void   RenderPalette(Float_t dx, Float_t dy) const;
protected:
  ITSModuleStepper* fM; // fModel dynamic-casted to ITSModuleStepperGL

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  ITSModuleStepperGL();
  virtual ~ITSModuleStepperGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  virtual Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(ITSModuleStepperGL, 0);
}; // endclass ITSModuleStepperGL

}

#endif
