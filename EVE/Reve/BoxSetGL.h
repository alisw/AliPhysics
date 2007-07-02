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

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  BoxSetGL();
  virtual ~BoxSetGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  virtual void Render(TGLRnrCtx & rnrCtx) { DirectDraw(rnrCtx); }

  ClassDef(BoxSetGL, 0);
}; // endclass BoxSetGL

}

#endif
