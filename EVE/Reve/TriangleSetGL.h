// $Header$

#ifndef REVE_TriangleSetGL_H
#define REVE_TriangleSetGL_H

#include <TGLObject.h>

class TGLRnrCtx;

namespace Reve {

class TriangleSet;

class TriangleSetGL : public TGLObject
{
private:
  TriangleSetGL(const TriangleSetGL&);            // Not implemented
  TriangleSetGL& operator=(const TriangleSetGL&); // Not implemented

protected:
  TriangleSet* fM; // Model object.

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  TriangleSetGL();
  virtual ~TriangleSetGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(TriangleSetGL, 0); // GL-renderer for TriangleSet class.
}; // endclass TriangleSetGL

}

#endif
