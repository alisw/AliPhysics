// $Header$

#ifndef REVE_BoxSetGL_H
#define REVE_BoxSetGL_H

#include <TGLObject.h>

namespace Reve {

class BoxSetGL : public TGLObject
{
protected:
  BoxSet* fM;

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  BoxSetGL();
  virtual ~BoxSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  virtual void Render(const TGLDrawFlags & flags) { DirectDraw(flags); }

  ClassDef(BoxSetGL, 0);
}; // endclass BoxSetGL

}

#endif
