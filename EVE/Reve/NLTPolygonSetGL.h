#ifndef REVE_NLTPolygonSetGL
#define REVE_NLTPolygonSetGL

#include <TGLObject.h>

namespace Reve {

class NLTPolygonSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  NLTPolygonSetGL();
  virtual  ~NLTPolygonSetGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }
  virtual void   SetBBox();

  ClassDef(NLTPolygonSetGL,0)  // GL renderer for NLTPolygonSet3D
};

}
#endif
