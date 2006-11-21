// $Header$

#ifndef REVE_QuadSetGL_H
#define REVE_QuadSetGL_H

#include <TGLObject.h>
#include <Reve/QuadSet.h>

namespace Reve {

class OldQuadSetGL : public TGLObject
{
protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

public:
  OldQuadSetGL();
  virtual ~OldQuadSetGL();

  virtual Bool_t SetModel(TObject* obj);
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

  virtual void DirectDraw(const TGLDrawFlags & flags) const;

  Bool_t SetupColor(const QuadSet::QuadBase& q) const;

  void   RenderQuads(const TGLDrawFlags & flags) const;
  void   RenderLines(const TGLDrawFlags & flags) const;

public:
  QuadSetGL();
  virtual ~QuadSetGL();

  virtual Bool_t SetModel(TObject* obj);
  virtual void   SetBBox();

  virtual Bool_t IgnoreSizeForOfInterest() const { return kTRUE; }

  ClassDef(QuadSetGL, 0);
};

}

#endif
