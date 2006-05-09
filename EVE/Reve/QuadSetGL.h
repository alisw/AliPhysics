// $Header$

#ifndef REVE_QuadSetGL_H
#define REVE_QuadSetGL_H

#include <TGLObject.h>

namespace Reve {

class QuadSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  QuadSetGL();
  virtual ~QuadSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  ClassDef(QuadSetGL, 0);
};

}

#endif
