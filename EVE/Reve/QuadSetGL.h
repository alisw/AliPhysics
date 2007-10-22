// $Header$

#ifndef REVE_QuadSetGL_H
#define REVE_QuadSetGL_H

#include <TGLObject.h>
#include <Reve/QuadSet.h>

namespace Reve {

class OldQuadSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

public:
  OldQuadSetGL();
  virtual ~OldQuadSetGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  ClassDef(OldQuadSetGL, 0);
};

/**************************************************************************/
/**************************************************************************/

class QuadSetGL : public TGLObject
{
  QuadSetGL(const QuadSetGL&);            // Not implemented
  QuadSetGL& operator=(const QuadSetGL&); // Not implemented

protected:
  QuadSet* fM;

  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

  Bool_t SetupColor(const DigitSet::DigitBase& q) const;

  void   RenderQuads(TGLRnrCtx & rnrCtx) const;
  void   RenderLines(TGLRnrCtx & rnrCtx) const;
  void   RenderHexagons(TGLRnrCtx & rnrCtx) const;

public:
  QuadSetGL();
  virtual ~QuadSetGL();

  virtual Bool_t ShouldDLCache(const TGLRnrCtx & rnrCtx) const;

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();

  virtual Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }

  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }
  virtual void   ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  ClassDef(QuadSetGL, 0);
};

}

#endif
