// $Header$

#ifndef REVE_LineGL_H
#define REVE_LineGL_H

#include <TGLObject.h>
#include <TPointSet3DGL.h>

class TGLViewer;
class TGLScene;

namespace Reve {

class Line;

class LineGL : public TPointSet3DGL
{
private:
  LineGL(const LineGL&);            // Not implemented
  LineGL& operator=(const LineGL&); // Not implemented

protected:
  Line* fM; // fModel dynamic-casted to LineGL

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  LineGL();
  virtual ~LineGL();

  virtual Bool_t SetModel(TObject* obj);

  // To support two-level selection
  // virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  // virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

  ClassDef(LineGL, 0);
}; // endclass LineGL

}

#endif
