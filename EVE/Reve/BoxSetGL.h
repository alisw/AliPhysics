// $Header$

#ifndef REVE_BoxSetGL_H
#define REVE_BoxSetGL_H

#include <TGLObject.h>
#include <Reve/BoxSet.h>

namespace Reve {

class BoxSetGL : public TGLObject
{
  BoxSetGL(const BoxSetGL&);            // Not implemented
  BoxSetGL& operator=(const BoxSetGL&); // Not implemented

protected:
  BoxSet* fM;

  mutable UInt_t  fBoxDL;

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

  Int_t  PrimitiveType() const;
  Bool_t SetupColor(const DigitSet::DigitBase& q) const;
  void   MakeOriginBox(Float_t* p, Float_t dx, Float_t dy, Float_t dz) const;
  void   RenderBox(const Float_t* p) const;
  void   MakeDisplayList() const;

public:
  BoxSetGL();
  virtual ~BoxSetGL();

  virtual Bool_t ShouldDLCache(const TGLRnrCtx & rnrCtx) const;
  virtual void   DLCacheDrop();
  virtual void   DLCachePurge();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void   ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  virtual void Render(TGLRnrCtx & rnrCtx);

  ClassDef(BoxSetGL, 0);
}; // endclass BoxSetGL

}

#endif
