#ifndef REVE_NLTPolygonSetGL
#define REVE_NLTPolygonSetGL

#ifndef ROOT_TGLObject
#include "TGLObject.h"
#endif


namespace Reve {

class NLTPolygonSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  NLTPolygonSetGL();
  virtual  ~NLTPolygonSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }
  virtual void   SetBBox();

  ClassDef(NLTPolygonSetGL,0)  // GL renderer for NLTPolygonSet3D
};

}
#endif
