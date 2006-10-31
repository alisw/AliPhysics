// $Header$

#ifndef REVE_QuadSetGL_H
#define REVE_QuadSetGL_H

#include <TGLObject.h>

namespace Reve {

class OldQuadSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  OldQuadSetGL();
  virtual ~OldQuadSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  ClassDef(OldQuadSetGL, 0);
};

}

#endif
