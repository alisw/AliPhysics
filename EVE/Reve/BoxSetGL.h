// $Header$

#ifndef REVE_BoxSetGL_H
#define REVE_BoxSetGL_H

#include <TGLObject.h>

namespace Reve {

class BoxSet;

class BoxSetGL : public TGLObject
{
  BoxSetGL(const BoxSetGL&);            // Not implemented
  BoxSetGL& operator=(const BoxSetGL&); // Not implemented

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
