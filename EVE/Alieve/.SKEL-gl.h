// $Header$

#ifndef ALIEVE_CLASS_H
#define ALIEVE_CLASS_H

#include <TGLObject.h>

namespace Alieve {

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
