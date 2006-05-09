// $Header$

#ifndef REVE_CLASS_H
#define REVE_CLASS_H

#include <TGLObject.h>

namespace Reve {

class CLASS
{
protected:
  STEM* fM; // fModel dynamic-casted to CLASS

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  CLASS();
  virtual ~CLASS();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  ClassDef(CLASS, 0);
}; // endclass CLASS

}

#endif
