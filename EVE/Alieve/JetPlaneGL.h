// $Header$

#ifndef ALIEVE_JetPlaneGL_H
#define ALIEVE_JetPlaneGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Alieve {

class JetPlane;

class JetPlaneGL : public TGLObject
{
private:
  JetPlaneGL(const JetPlaneGL&);            // Not implemented
  JetPlaneGL& operator=(const JetPlaneGL&); // Not implemented

protected:
  JetPlane* fM; // fModel dynamic-casted to JetPlaneGL

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  JetPlaneGL();
  virtual ~JetPlaneGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  // To support two-level selection
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(JetPlaneGL, 0);
}; // endclass JetPlaneGL

}

#endif
