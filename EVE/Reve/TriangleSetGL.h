// $Header$

#ifndef REVE_TriangleSetGL_H
#define REVE_TriangleSetGL_H

#include <TGLObject.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class TriangleSet;

class TriangleSetGL : public TGLObject
{
private:
  TriangleSetGL(const TriangleSetGL&);            // Not implemented
  TriangleSetGL& operator=(const TriangleSetGL&); // Not implemented

protected:
  TriangleSet* fM; // fModel dynamic-casted to TriangleSetGL

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  TriangleSetGL();
  virtual ~TriangleSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(TriangleSetGL, 0);
}; // endclass TriangleSetGL

}

#endif
